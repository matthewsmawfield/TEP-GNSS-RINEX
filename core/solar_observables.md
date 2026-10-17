# Coherent two-way and three-way observations

`solar_observables.py` predicts round-trip light time and count-averaged coherent Doppler. A two-way observation uses one ground station for transmission and reception. A three-way observation uses separate stations. The spacecraft trajectory, ground station trajectories and clock transformations are supplied explicitly by callbacks.

The Doppler convention is receiver reference phase minus received coherent phase, divided by the reception count interval:

\[
F_D=\frac{1}{T_c}\left[\int_{\tau_{R,a}}^{\tau_{R,b}}f_{\rm ref}(\tau_R)d\tau_R
-M_2\int_{\tau_{T,a}}^{\tau_{T,b}}f_T(\tau_T)d\tau_T\right].
\]

The retarded downlink and uplink solutions determine the transmission interval. Separate station clocks distinguish two-way and three-way measurements. Transmitted frequency ramps, receiver reference ramps and constant receiver references are supported. Integration splits at every relevant ramp boundary. Missing coverage or conflicting overlapping ramps fail explicitly.

The implementation evaluates the phase-rate difference directly rather than subtracting large accumulated phases or nearby floating-point timestamps. Frequency arithmetic uses Decimal, and radio turnaround ratios use exact fractions. Cassini Ka/Ka uses 14/15; X/X uses 880/749 and X/Ka uses 3344/749.

`CoherentModel.round_trip` returns both retarded epochs and coordinate round-trip delay, including explicit antenna/electronics and transponder delays. `CoherentModel.doppler` returns the predicted Doppler, transmission count width and endpoint delay diagnostics. The Sun and other bodies can contribute first-PN Shapiro delay. A delay callback can add a chosen propagation model separately to each leg.

For a TEP calculation, a propagation callback supplies the ray contribution of the selected field solution. TEP changes to endpoint clock response and spacecraft motion must also be supplied through the clock and state callbacks. A static ray contribution alone does not supply those effects.

The SCE1 adapter and archive checks are in `notes/solar-observation-2026-10-12/`. They use label-decoded ODF records, station-specific frequency ramps and SPICE ephemerides. Archived-versus-modern Earth orientation/ephemeris conventions, tropospheric baseline and increments, ionosphere and solar plasma, antenna/transponder calibration, and trajectory/nuisance fitting still require reconciliation before a precision gravity or TEP likelihood can be reported. The current real-data pilot is an uncalibrated operator check.

Run the analytic and precision tests with:

```sh
python TEP/tests/test_solar_observables.py -v
```

The observation convention follows [Moyer, *Formulation for Observed and Computed Values of Deep Space Network Data Types for Navigation*, equation 13–50](https://descanso.jpl.nasa.gov/monograph/series2/Descanso2_all.pdf). Data are from the [NASA/PDS Cassini SCE1 archive](https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=CO-SS-RSS-1-SCE1-V1.0).
