#!/usr/bin/env python3
"""Background solver for the attractorE / attractorT CosmoLattice models.

Given (model, alphaAtt, k, kappa, N_*), integrates the homogeneous inflaton
equation in e-folds, finds the end of inflation (epsilon_H = 1) and the field
value N_* e-folds earlier, fixes lambda from A_s, and prints (or writes) the
CosmoLattice parameters: lambda and the initial amplitude/momentum at the end
of inflation in GeV / GeV^2.

Potentials (x = phi/MPl, b = sqrt(2/(3 alpha)), c = b/2), V = (3/4) lambda MPl^4 F(x):
  E: F = [kappa (1 - cosh bx) + sinh bx]^k
  T: F = [(1 + kappa + (1 - kappa) cosh bx)/2]^2 tanh^k(cx)
See notes/inflation_models.tex.
"""
import argparse
import numpy as np
from scipy.integrate import solve_ivp

MPL = 2.435e18  # reduced Planck mass in GeV
AS = 2.1e-9


def shape(model, alpha, k, kappa):
    b = np.sqrt(2 / (3 * alpha))
    c = b / 2
    if model == "E":
        F = lambda x: (np.sinh(b * x) - 2 * kappa * np.sinh(b * x / 2) ** 2) ** k
        dF = lambda x: k * (np.sinh(b * x) - 2 * kappa * np.sinh(b * x / 2) ** 2) ** (k - 1) \
            * b * (np.cosh(b * x) - kappa * np.sinh(b * x))
    else:
        D = lambda x: 1 + (1 - kappa) * np.sinh(b * x / 2) ** 2
        dD = lambda x: (1 - kappa) * b * np.sinh(b * x) / 2
        F = lambda x: D(x) ** 2 * np.tanh(c * x) ** k
        dF = lambda x: 2 * D(x) * dD(x) * np.tanh(c * x) ** k \
            + D(x) ** 2 * k * c * np.tanh(c * x) ** (k - 1) / np.cosh(c * x) ** 2
    return F, dF


def solve(model, alpha, k, kappa, nstar):
    F, dF = shape(model, alpha, k, kappa)

    # Start deep on the plateau, on the slow-roll attractor, and integrate in N.
    # y = (x, dx/dN); x'' = -(3 - eps)(x' + F'/F), eps = x'^2/2.
    def rhs(N, y):
        x, p = y
        eps = p * p / 2
        return [p, -(3 - eps) * (p + dF(x) / F(x))]

    def end(N, y):
        return y[1] ** 2 / 2 - 1
    end.terminal = True
    end.direction = 1

    # Pick x0 giving comfortably more than nstar + 10 e-folds (slow-roll estimate).
    x0 = 1.0
    while True:
        xs = np.linspace(0.05, x0, 4000)
        nsr = np.trapezoid(F(xs) / dF(xs), xs)
        if nsr > nstar + 15:
            break
        x0 *= 1.2
    y0 = [x0, -dF(x0) / F(x0)]
    sol = solve_ivp(rhs, [0, 1e4], y0, events=end, dense_output=True, rtol=1e-10, atol=1e-12)
    Nend = sol.t_events[0][0]
    xend, pend = sol.y_events[0][0]

    Nst = Nend - nstar
    xst, pst = sol.sol(Nst)
    epsV = 0.5 * (dF(xst) / F(xst)) ** 2
    # A_s = V/(24 pi^2 eps_V MPl^4) with V = (3/4) lambda MPl^4 F
    lam = 32 * np.pi ** 2 * AS * epsV / F(xst)

    # dphi/dt = MPl H x',  H^2 = V/(MPl^2 (3 - eps_H))
    H = MPL * np.sqrt(0.75 * lam * F(xend) / (3 - pend ** 2 / 2))
    phidot = MPL * H * pend

    # Slow-roll observables at x_*
    d2 = lambda x, h=1e-4: (F(x + h) - 2 * F(x) + F(x - h)) / h ** 2
    etaV = d2(xst) / F(xst)
    ns = 1 - 6 * epsV + 2 * etaV
    r = 16 * epsV
    return dict(lam=lam, xend=xend, xstar=xst, phi0=xend * MPL, pi0=phidot, ns=ns, r=r, H=H)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", choices=["E", "T"])
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--k", type=int, default=2)
    ap.add_argument("--kappa", type=float, default=1.0)
    ap.add_argument("--nstar", type=float, default=55.0)
    args = ap.parse_args()
    if args.k < 2 or args.k % 2:
        ap.error("k must be an even integer >= 2")
    s = solve(args.model, args.alpha, args.k, args.kappa, args.nstar)
    print(f"# {args.model}-model alphaAtt={args.alpha} k={args.k} kappa={args.kappa} N*={args.nstar}")
    print(f"# slow roll at x_*={s['xstar']:.4f}: ns={s['ns']:.4f} r={s['r']:.5f};  x_end={s['xend']:.4f}, H_end={s['H']:.4e} GeV")
    print(f"initial_amplitudes = {s['phi0']:.6e} 0 # GeV")
    print(f"initial_momenta = {s['pi0']:.6e} 0 # GeV^2")
    print(f"lambda = {s['lam']:.6e}")
    print(f"alphaAtt = {args.alpha}")
    print(f"k = {args.k}")
    print(f"kappa = {args.kappa}")


if __name__ == "__main__":
    main()
