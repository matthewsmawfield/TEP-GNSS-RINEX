"""Coherent two-/three-way radiometric operator.

Station-time phase counting follows Moyer eq. 13-50. Ephemeris callbacks
return barycentric km and km/s at TDB seconds past J2000. Clock callbacks
return station-time minus TDB seconds. Media/hardware/TEP delays enter as
explicit leg callbacks; an omitted callback is not a certified correction.
"""
from dataclasses import dataclass
from typing import Callable
import numpy as np
from decimal import Decimal,localcontext
from fractions import Fraction
from scipy.optimize import brentq

C_KM_S=299792.458
LD=np.longdouble
def dec(value):
    if isinstance(value,Decimal):return value
    if isinstance(value,Fraction):return Decimal(value.numerator)/Decimal(value.denominator)
    return Decimal(str(value))

@dataclass(frozen=True)
class Ramp:
    start: float
    end: float
    frequency_hz: float
    rate_hz_s: float

class RampTable:
    def __init__(self,ramps):
        self.ramps=sorted(ramps,key=lambda r:r.start)
        if any(r.end<=r.start for r in self.ramps):raise ValueError('Invalid ramp interval')
        for i,a in enumerate(self.ramps):
            for b in self.ramps[i+1:]:
                if b.start>=a.end:break
                for time in [max(a.start,b.start),min(a.end,b.end)]:
                    fa=dec(a.frequency_hz)+dec(a.rate_hz_s)*(dec(time)-dec(a.start))
                    fb=dec(b.frequency_hz)+dec(b.rate_hz_s)*(dec(time)-dec(b.start))
                    if abs(fa-fb)>Decimal('0.000001'):
                        raise ValueError('Conflicting overlapping frequency ramps')
    def frequency(self,time):
        t=dec(time)
        for r in self.ramps:
            if dec(r.start)<=t<=dec(r.end):
                return dec(r.frequency_hz)+dec(r.rate_hz_s)*(t-dec(r.start))
        raise ValueError('Transmit ramp missing at emission epoch')
    def integral(self,start,width):
        """Integrate locally, never subtract huge accumulated phases."""
        if width<=0:raise ValueError('Positive count width required')
        with localcontext() as context:
            context.prec=50
            cursor=dec(start);end=cursor+dec(width);cycles=Decimal(0)
            for r in self.ramps:
                a,b=dec(r.start),dec(r.end)
                if b<=cursor:continue
                if a>cursor:raise ValueError('Uncovered transmit ramp interval')
                stop=min(end,b);dt=stop-cursor
                f=dec(r.frequency_hz)+dec(r.rate_hz_s)*(cursor-a)
                cycles+=f*dt+dec(r.rate_hz_s)*dt*dt/2
                cursor=stop
                if cursor==end:return cycles
            raise ValueError('Transmit ramp coverage ends before count interval')

def cassini_turnaround(up,down):
    if up==3 and down==3:return Fraction(14,15)
    numerators={1:240,2:880,3:3344}
    denominators={1:221,2:749,3:3599}
    if up not in denominators or down not in numerators:raise ValueError('Unsupported radio band')
    return Fraction(numerators[down],denominators[up])

def reference_multiplier(exciter,down):
    # ODF constant reference is represented in the declared exciter band.
    # Ka/Ka Cassini uses its nonstandard ratio; S/X references use normal
    # ratios. Ramped receiver-reference data need a reference phase callback.
    return cassini_turnaround(exciter,down)

def shapiro_delay(emitter,receiver,body_position,gm,gamma=1.):
    a=np.asarray(emitter,dtype=LD)-np.asarray(body_position,dtype=LD)
    b=np.asarray(receiver,dtype=LD)-np.asarray(body_position,dtype=LD)
    ra=np.sqrt(np.dot(a,a));rb=np.sqrt(np.dot(b,b))
    vector=b-a;distance=np.sqrt(np.dot(vector,vector))
    denominator=ra+rb-distance
    if denominator<=0:raise ValueError('Singular or occulted Shapiro geometry')
    return LD(1+gamma)*LD(gm)/LD(C_KM_S)**3*np.log1p(2*distance/denominator)

class CoherentModel:
    def __init__(self,state:Callable,clock=None,bodies=(),gamma=1.,extra_delay=None,clock_rate=None):
        self.state=state
        self.clock=clock or (lambda station,t:0.)
        self.bodies=list(bodies)
        self.gamma=gamma
        self.extra_delay=extra_delay
        self.clock_rate=clock_rate or (lambda station,t:(self.clock(station,float(t)+2)-self.clock(station,float(t)-2))/4)
        self.quadrature_order=5
    def to_coordinate(self,station,station_time):
        t=float(station_time)
        for _ in range(12):
            new=float(LD(station_time)-LD(self.clock(station,t)))
            if abs(new-t)<1e-8:return LD(station_time)-LD(self.clock(station,new))
            t=new
        raise RuntimeError('Station-time conversion did not converge')
    def position(self,target,t):
        """Apply an epoch residual when the platform retains extended precision.

        Doppler precision relies on phase-rate arithmetic, including on
        platforms where longdouble has the same precision as float64.
        """
        base=float(t)
        state=np.asarray(self.state(target,base),dtype=LD)
        return state[:3]+state[3:6]*(LD(t)-LD(base))
    def leg_rate_delta(self,emitter,receiver,emission_t,receive_t,leg_name):
        """d(emission time)/d(reception time)-1, avoiding near-unit subtraction."""
        pe=self.position(emitter,emission_t);pr=self.position(receiver,receive_t)
        ve=np.asarray(self.state(emitter,float(emission_t))[3:6],dtype=float)
        vr=np.asarray(self.state(receiver,float(receive_t))[3:6],dtype=float)
        ray=np.asarray(pr-pe,dtype=float);distance=np.linalg.norm(ray);n=ray/distance
        de,dr=0.,0.
        middle=float((emission_t+receive_t)/2)
        for body,gm in self.bodies:
            pb=np.asarray(self.position(body,middle),dtype=float)
            vb=np.asarray(self.state(body,middle)[3:6],dtype=float)
            a=np.asarray(pe,dtype=float)-pb;b=np.asarray(pr,dtype=float)-pb
            ra=np.linalg.norm(a);rb=np.linalg.norm(b);s=ra+rb
            dse=np.dot(a/ra,ve-vb/2)+np.dot(b/rb,-vb/2)
            dsr=np.dot(a/ra,-vb/2)+np.dot(b/rb,vr-vb/2)
            dre=-np.dot(n,ve);drr=np.dot(n,vr)
            coeff=(1+self.gamma)*gm/C_KM_S**3
            de+=coeff*((dse+dre)/(s+distance)-(dse-dre)/(s-distance))
            dr+=coeff*((dsr+drr)/(s+distance)-(dsr-drr)/(s-distance))
        if self.extra_delay:
            h=.2
            def correction(te,tr):
                return float(self.extra_delay(emitter,receiver,te,tr,self.position(emitter,te),self.position(receiver,tr),leg_name))
            de+=(correction(emission_t+h,receive_t)-correction(emission_t-h,receive_t))/(2*h)
            dr+=(correction(emission_t,receive_t+h)-correction(emission_t,receive_t-h))/(2*h)
        a=float(np.dot(n,vr)/C_KM_S)+dr
        b=float(np.dot(n,ve)/C_KM_S)-de
        return (b-a)/(1-b)
    def leg(self,emitter,receiver,receive_t,leg_name):
        pr=self.position(receiver,receive_t)
        initial=self.position(emitter,receive_t)
        delay=np.linalg.norm(pr-initial)/LD(C_KM_S)
        for _ in range(30):
            emission_t=LD(receive_t)-delay
            pe=self.position(emitter,emission_t)
            new=np.linalg.norm(pr-pe)/LD(C_KM_S)
            middle=(emission_t+receive_t)/2
            for body,gm in self.bodies:
                pb=self.position(body,middle)
                new+=shapiro_delay(pe,pr,pb,gm,self.gamma)
            if self.extra_delay:
                new+=LD(self.extra_delay(emitter,receiver,emission_t,receive_t,pe,pr,leg_name))
            # On platforms where longdouble is float64, an AU-scale
            # delay has representable steps larger than 1e-13 seconds.
            # Include the propagated emitter-epoch quantisation error,
            # so a round-off two-cycle is not called non-convergence.
            velocity=np.linalg.norm(np.asarray(self.state(emitter,float(emission_t))[3:6],dtype=LD))
            tolerance=max(LD('1e-13'),4*abs(np.spacing(LD(new))),
                          2*abs(np.spacing(LD(emission_t)))*velocity/LD(C_KM_S))
            if abs(new-delay)<=tolerance:
                return LD(receive_t)-new,new
            delay=new
        raise RuntimeError('Retarded light-time iteration did not converge')
    def round_trip(self,tx,rx,spacecraft,receive_t,transponder_delay_s=0.,
                   uplink_delay_s=0.,downlink_delay_s=0.):
        antenna_receive=receive_t-downlink_delay_s
        bounce,down=self.leg(spacecraft,rx,antenna_receive,'downlink')
        uplink_receive=bounce-transponder_delay_s
        emit,up=self.leg(tx,spacecraft,uplink_receive,'uplink')
        return {'emit_t':emit-uplink_delay_s,'bounce_t':bounce,'receive_t':receive_t,
                'coordinate_delay_s':up+down+LD(transponder_delay_s)+LD(uplink_delay_s)+LD(downlink_delay_s)}
    def doppler(self,tx,rx,spacecraft,midpoint_station_s,count_s,ramps,
                turnaround,reference_frequency_hz,reference_ratio=1.,
                transponder_delay_s=0.,reference_phase=None,
                uplink_delay_s=0.,downlink_delay_s=0.,reference_ramps=None):
        """Moyer ramped F2/F3: receiver reference minus received phase/TC.

        Separate stations distinguish F3 from F2. Clock and precision delay
        differences determine the transmit width, rather than subtraction of
        two large floating-point emission epochs.
        """
        if count_s<=0:raise ValueError('Invalid count interval')
        if reference_phase is not None and reference_ramps is not None:
            raise ValueError('Choose one receiver reference phase representation')
        sa=LD(midpoint_station_s)-LD(count_s)/2
        sb=LD(midpoint_station_s)+LD(count_s)/2
        ta=self.to_coordinate(rx,sa);tb=self.to_coordinate(rx,sb)
        a=self.round_trip(tx,rx,spacecraft,ta,transponder_delay_s,uplink_delay_s,downlink_delay_s)
        b=self.round_trip(tx,rx,spacecraft,tb,transponder_delay_s,uplink_delay_s,downlink_delay_s)
        txa=a['emit_t'];txb=b['emit_t']
        d_rx=LD(self.clock(rx,tb))-LD(self.clock(rx,ta))
        d_tx=LD(self.clock(tx,txb))-LD(self.clock(tx,txa))
        width=LD(count_s)-(b['coordinate_delay_s']-a['coordinate_delay_s'])+d_tx-d_rx
        start=LD(ta)-a['coordinate_delay_s']+LD(self.clock(tx,txa))
        # Keep endpoint coverage checks, but integrate phase rate to avoid
        # finite differencing kilometre-scale delays at a ~1e8-second epoch.
        # Decimal emission endpoints avoid a false gap at a ramp knot.
        # Full coverage is checked below after their station-time conversion.
        # Split receiver integration where an emission interval crosses a
        # ramp boundary. One quadrature across a rate/frequency break is wrong.
        def emission_station(st):
            tr=self.to_coordinate(rx,st)
            trip=self.round_trip(tx,rx,spacecraft,tr,transponder_delay_s,uplink_delay_s,downlink_delay_s)
            return dec(st)-dec(trip['coordinate_delay_s'])+dec(self.clock(tx,trip['emit_t']))-dec(self.clock(rx,tr))
        sa_dec=dec(midpoint_station_s)-dec(count_s)/2
        ea=emission_station(sa_dec);eb=emission_station(sa_dec+dec(count_s))
        ramps.integral(ea,eb-ea)
        boundaries=sorted({dec(t) for r in ramps.ramps for t in [r.start,r.end] if ea<dec(t)<eb})
        cuts=[0.]
        for boundary in boundaries:
            root=brentq(lambda dt:float(emission_station(sa_dec+dec(dt))-boundary),0.,float(count_s),xtol=1e-9)
            cuts.append(root)
        if reference_ramps is not None:
            reference_ramps.integral(sa_dec,dec(count_s))
            cuts.extend(float(dec(t)-sa_dec) for r in reference_ramps.ramps
                        for t in [r.start,r.end] if sa_dec<dec(t)<sa_dec+dec(count_s))
        cuts=sorted(set(cuts))
        cuts.append(float(count_s))
        base_nodes,base_weights=np.polynomial.legendre.leggauss(self.quadrature_order)
        node_times=[];weights=[]
        for lo,hi in zip(cuts,cuts[1:]):
            for node,w in zip(base_nodes,base_weights):
                node_times.append(sa_dec+dec((lo+hi)/2)+dec(node*(hi-lo)/2))
                weights.append(w*(hi-lo)/float(count_s))
        values=[];rates=[]
        with localcontext() as context:
            context.prec=50
            ref=dec(reference_ratio)*dec(reference_frequency_hz)
            for st in node_times:
                tr=self.to_coordinate(rx,st)
                trip=self.round_trip(tx,rx,spacecraft,tr,transponder_delay_s,uplink_delay_s,downlink_delay_s)
                tx_time=trip['emit_t'];bounce=trip['bounce_t']
                down=self.leg_rate_delta(spacecraft,rx,bounce,tr-downlink_delay_s,'downlink')
                up=self.leg_rate_delta(tx,spacecraft,tx_time+uplink_delay_s,bounce-transponder_delay_s,'uplink')
                geometric=down+up+down*up
                rtx=float(self.clock_rate(tx,tx_time));rrx=float(self.clock_rate(rx,tr))
                clocks=(rtx-rrx)/(1+rrx)
                rate_delta=geometric+clocks+geometric*clocks
                emit_st=st-dec(trip['coordinate_delay_s'])+dec(self.clock(tx,tx_time))-dec(self.clock(rx,tr))
                frequency=ramps.frequency(emit_st)
                received=dec(turnaround)*frequency
                reference=ref if reference_ramps is None else dec(reference_ratio)*reference_ramps.frequency(st)
                values.append(reference-received-received*dec(rate_delta))
                rates.append(rate_delta)
            if reference_phase is not None:
                extra=dec(reference_phase(sa_dec,sa_dec+dec(count_s)))/dec(count_s)-ref
                values=[v+extra for v in values]
            value=sum(dec(w)*v for w,v in zip(weights,values))/2
        stable_width=count_s*(1+np.dot(weights,rates)/2)
        return {'doppler_hz':float(value),'transmit_count_width_s':float(stable_width),
                'round_trip_start_s':float(a['coordinate_delay_s']),
                'round_trip_end_s':float(b['coordinate_delay_s']),
                'coordinate_receive_start':ta,'coordinate_receive_end':tb,
                'tx_start_station_s':float(start)}
