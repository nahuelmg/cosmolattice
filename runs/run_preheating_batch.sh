#!/usr/bin/env bash
# Preheating batch for the benchmark points of notes/inflation_models.tex ("Recommended minimal set of runs"),
# with the daughter coupling (1/2) g^2 phi^2 chi^2 switched on and GWs. Outputs in runs/benchmarks/<name>/.
# g = 2 sqrt(q0) omega_*/phi_i, with q0 = g^2 phi_i^2/(4 omega_*^2) the initial resonance parameter.
# q grows as a^((6k-24)/(k+2)): constant for k=4, ~a^1.5 for k=6, ~a^3 for k=10. q0 is therefore lowered
# for k>=6 so that chi's program-time frequency (~2 sqrt(q)) stays resolved by dt until the end of the run.
# For k=10 (alpha_t=2) the wave speed in program time grows like a, so dt is small and tMax stops near a~20-25.
# The undeformed k=2 points (Starobinsky E, T k=2) are runs/E_starobinsky_gw and runs/T_k2_a1_gw (same settings).
set -euo pipefail
R="$(dirname "${BASH_SOURCE[0]}")/run_benchmarks.sh"
COMMON=(kIR=0.5 kCutOff=20 withGWs=true tOutputInfreq=10)

# k = 4: q0 = 1e4 (q constant); dt = 0.005 because chi oscillates at ~2 sqrt(q) = 200 in program time all run long
# (dt = 0.01 sits at the leapfrog limit omega dt = 2: Friedmann violation grew to 4.5e-2, see benchmarks/E_k4_a1_dt0p01)
"$R" E_k4_a1              -- "${COMMON[@]}" dt=0.005 tMax=300 g=2.923e-03
"$R" E_k4_a1_kappa0p99999 -- "${COMMON[@]}" dt=0.005 tMax=300 g=3.150e-03
"$R" T_k4_a1_kappa0p99999 -- "${COMMON[@]}" dt=0.005 tMax=300 g=7.749e-04
# k = 6, 10: NOT RUN. The quadratic coupling does not give broad resonance for these benchmarks.
# In the E k=6 alpha=1 run with q0=10 (benchmarks/E_k6_a1_q10_kIR0p5) chi never grew: the real
# oscillation period is ~20 program-time units (not 2 pi), so q_eff is far below the estimate, and the
# only unstable modes are k <~ 0.5 at the lattice IR edge. A single-mode scan on that background gives
# n growth of only 1e3 / 6e3 / 1e5 by t=150 for g x10 / x30 / x100, with chi frequencies 800-8000 in
# program time (dt <~ 1e-3 - 1e-4). Not feasible here; the lines are kept for reference.
# k = 6: q0 = 10 (grows to ~3e3 by a~45)
# "$R" E_k6_a1              -- "${COMMON[@]}" dt=0.008 tMax=300 g=1.322e-04
# "$R" E_k6_a5              -- "${COMMON[@]}" dt=0.008 tMax=300 g=3.723e-05
# "$R" E_k6_a4_kappa0p9999  -- "${COMMON[@]}" dt=0.008 tMax=300 g=5.080e-05
# "$R" T_k6_a1              -- "${COMMON[@]}" dt=0.008 tMax=300 g=2.258e-05
# "$R" T_k6_a1_kappa0p99999 -- "${COMMON[@]}" dt=0.008 tMax=300 g=2.428e-05
# k = 10: q0 = 1 (grows ~a^3), Courant-limited
# "$R" E_k10_a5             -- "${COMMON[@]}" dt=0.0025 tMax=75 tOutputInfreq=5 g=2.627e-05
# "$R" T_k10_a1             -- "${COMMON[@]}" dt=0.0025 tMax=45 tOutputInfreq=3 g=1.299e-05
# deformed k = 2: q0 = 1e4
"$R" E_k2_a1_kappa0p9999  -- "${COMMON[@]}" dt=0.01 tMax=300 g=5.260e-03
"$R" T_k2_a1_kappa0p9999  -- "${COMMON[@]}" dt=0.01 tMax=300 g=1.987e-03
