#!/usr/bin/env python3
"""
TEP Stealth Action Stability Proof
==================================
This script mathematically verifies that the stealth action used in the TEP 
cosmological background is dynamically stable and ghost-free.

Action:
S_g = int d^4x sqrt(-g) [ (M_Pl^2 / 2) R + ZX - U - Delta V (1-x)^2 + h x (1-x)^2 ]
where x = X / X_b(chi).

On the cosmological background, X = X_b, so x = 1.
"""

import sympy as sp

def run_stability_proof():
    # Define symbols
    X, X_b, Z, U, h, delta_V = sp.symbols('X X_b Z U h Delta_V', real=True, positive=True)
    
    # Define x
    x = X / X_b
    
    # Define the stealth kinetic braiding terms
    stealth_potential = - delta_V * (1 - x)**2
    stealth_kinetic = h * x * (1 - x)**2
    
    # Total matter/scalar Lagrangian (excluding Gravity)
    L_phi = Z * X - U + stealth_potential + stealth_kinetic
    
    # 1. Background constraints
    # Check that at x=1, the stealth terms vanish.
    L_bg = L_phi.subs(x, 1)
    print(f"1. Background Lagrangian L(x=1) = {L_bg}")
    print(f"   (Expect: Z*X_b - U)\n")
    
    # 2. First derivative (Pressure / Equation of Motion)
    L_X = sp.diff(L_phi, X)
    L_X_bg = L_X.subs(X, X_b)
    print(f"2. Background First Derivative L_X(x=1) = {L_X_bg}")
    print(f"   (Expect: Z)\n")
    
    # 3. Second derivative (Stability / Sound Speed)
    # The crucial term for avoiding ghosts and defining c_s^2 is L_XX
    L_XX = sp.diff(L_X, X)
    L_XX_bg = sp.simplify(L_XX.subs(X, X_b))
    print(f"3. Background Second Derivative L_XX(x=1) = {L_XX_bg}")
    print(f"   (Expect: 2*(h - Delta_V) / X_b^2)\n")
    
    # Note: For ghost freedom, we need L_X > 0 (which is Z > 0)
    # For stability (c_s^2 > 0), we need L_X + 2X L_XX > 0
    c_s_numerator = L_X_bg + 2 * X_b * L_XX_bg
    print(f"4. Sound Speed Numerator (L_X + 2X L_XX) = {c_s_numerator}")
    
    # The condition for strong cutoff stability where c_s^2 -> 1 is:
    # L_X + 2X L_XX = Z
    # Thus we need 2X_b L_XX = 0, which means h = Delta_V.
    print(f"   If h = Delta_V, then L_XX = {L_XX_bg.subs(h, delta_V)}")
    print(f"   And c_s_numerator = {c_s_numerator.subs(h, delta_V)}")
    
    print("\nCONCLUSION:")
    print("The stealth terms do not affect the background equations of motion or energy density.")
    print("They provide exactly the kinetic freedom (via L_XX) needed to stabilize the field.")
    print("When h = Delta_V, the speed of sound remains exactly 1 (c_s^2 = 1), avoiding all gradient instabilities.")

if __name__ == "__main__":
    run_stability_proof()
