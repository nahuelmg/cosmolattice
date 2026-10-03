# Summary line for finished runs: python3 summarize.py <run dir> ...
import sys, numpy as np
for d in sys.argv[1:]:
    try:
        A = np.loadtxt(d + "/average_energies.txt"); t = A[:, 0]
        K0, G0, K1, G1, V0, Vi, Et = [A[:, i] for i in range(1, 8)]
        w = (K0 + K1 - (G0 + G1) / 3 - V0 - Vi) / Et
        m = t > t[-1] - 20
        c = np.loadtxt(d + "/average_energy_conservation.txt")
        a = np.loadtxt(d + "/average_scale_factor.txt")
        g = np.loadtxt(d + "/average_energies_gws.txt")
        print(f"{d}: t={t[-1]:.0f} a={a[-1,1]:.1f} <w>_last20={w[m].mean():.3f} "
              f"chi_frac={((K1+G1)/Et)[m].mean():.2f} Vint={(Vi/Et)[m].mean():.2f} "
              f"cons_max={np.abs(c[:,1]).max():.1e} rhoGW/rho={g[-1,1]:.1e} nan={bool(np.isnan(A).any())}")
    except Exception as e:
        print(f"{d}: summary failed: {e}")
