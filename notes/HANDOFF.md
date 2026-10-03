# Handoff: α-attractor models in CosmoLattice

Read this first. It is the state at the end of the session of 2026-10-02.

**This repo is used on several computers.** All paths below are relative to the repository root
(the directory containing `CMakeLists.txt`); never assume a fixed home directory, user name or conda location.
Check `git remote -v` and the available tools on each machine before using them.

## Goal of the work
Implement in CosmoLattice the inflaton potentials that remain viable after the latest CMB data, following
**arXiv:2510.18656** (Ellis, Garcia, Olive, Verner: attractor models vs Planck/BK18/ACT DR6/SPT-3G).
Then simulate preheating with GWs and redshift the GW spectra to today, following **arXiv:1707.04533**
(Figueroa & Torrenti). Both papers' LaTeX sources are in `papers/`.

## Repository state
- Repo: CosmoLattice **2.0** (upstream commit `acc8278d`), cloned at a different place on each computer. In 2.0, models live in `models/`, not `src/models/`.
- Work branch: **`attractor-models`**, pushed to the user's fork `git@github.com:nahuelmg/cosmolattice.git` (SSH works; HTTPS has no saved credentials).
- **Remote names differ between computers.** On one machine the fork is `fork` and `origin` is the official `cosmolattice/cosmolattice`; on another, `origin` is the fork itself. Run `git remote -v` and push only to the remote whose URL is `nahuelmg/cosmolattice`. **Never push to `cosmolattice/cosmolattice`.** Local `master` = upstream master, untouched.
- Machine-local state (not in git, may not exist on the current computer): `git stash@{0}` ("local changes before pulling v2.0") with the user's old edits to the v1 files `src/models/parameter-files/lphi4.in` (deleted) and `lphi4U1.in`; the `build_attractor{E,T}/` directories.
- **Missing files:** on the first computer, the user's untracked v1 files `src/models/lphi4_only.h`, `src/models/m2phi2.h` and `src/models/parameter-files/lphi4_{1..5,4bis,only}.in`, `m2phi2.in` disappeared from disk mid-session (around 18:41). The cause is unknown; they are not in git or the Trash. The user was told. If they recover them, they may want them ported to v2.0.
- `gh` CLI is not installed: no PR creation from here. PR link: https://github.com/nahuelmg/cosmolattice/pull/new/attractor-models

## What exists (all on branch `attractor-models`)
| Path | What |
|---|---|
| `models/attractorE.h` | E-model, generalized + deformed: `V = (3/4)λ MP⁴ [κ(1−cosh bx)+sinh bx]^k`, x=φ/MP, b=√(2/3α) |
| `models/attractorT.h` | T-model, generalized + deformed: `V = (3/4)λ MP⁴ D² tanh^k(cx)`, D=1+(1−κ)sinh²(bx/2), c=b/2 |
| `models/parameter-files/attractor{E,T}_*.in` | 14 benchmark points of the paper (λ and initial conditions from the helper) |
| `notes/attractor_ics.py` | Background solver: `python attractor_ics.py E --alpha 1 --k 2 --kappa 1 --nstar 50`, which prints λ (from A_s) and `initial_amplitudes`/`initial_momenta` at ε_H=1 |
| `notes/inflation_models.tex/.pdf` | Physics note: data, viable models, V, V′, V″, benchmarks table |
| `notes/implementation_report.tex/.pdf` | How it was built, files, parameter choices, validation, and the derivation of the GW redshift to today (f, h²Ω_GW) |
| `runs/E_starobinsky_gw/`, `runs/T_k2_a1_gw/` | Preheating runs with GWs (`run.in` + full output) |
| `runs/attractor_runs.ipynb` | Executed notebook: energies, w, field/GW spectra, Sec. 5 with today's f and h²Ω_GW, Sec. 6 the benchmark set |
| `runs/run_benchmarks.sh` | Generic driver: copies `models/parameter-files/attractor<name>.in` to `runs/benchmarks/<name>/run.in`, applies `-- key=value` overrides, runs, skips finished runs (`done` file) |
| `runs/run_preheating_batch.sh` | The benchmark preheating batch (per-run g, dt, tMax, with the reasoning in comments) |
| `runs/benchmarks/` | Outputs of that batch + reference/failed runs (`E_k4_a1_dt0p01`, `E_k6_a1_q10_kIR0p5`, `selfres_E_k6_a1`), `batch.log`, `summarize.py` |
| `runs/make_notebook.py` | Generates the notebook. **Edit this, not the .ipynb**, then rebuild (see below) |

### Model header conventions
- Run-time parameters: `lambda`, `alphaAtt` (NOT `alpha`, which is CosmoLattice's time-rescaling exponent), `k` (even integer ≥2, checked), `kappa` (1 = undeformed), `g` (daughter coupling ½g²φ²χ², default 0).
- 2 scalars: inflaton (0) and daughter χ (1). The potential is in program units via `pref = ¾λ MP²/ω*²`.
- Program variables: `fStar = MPl`, `omegaStar = sqrt(λ_k) MPl x_i^{(k−2)/2}` with `λ_k = k·¾λ·b^k` (E) or `c^k` (T), and `alpha = 3(k−2)/(k+2)`. For k=2: ω* = inflaton mass, α=0. For k=4: α=1.
- `initial_amplitudes` are in GeV, `initial_momenta` in GeV².

## How to build / run / rebuild the notebook
Tools: needs cmake, a C++ compiler, FFTW3, and Python 3 with numpy/scipy/matplotlib/nbformat/jupyter.
Where they come from depends on the computer: on the first one, a conda env `cosmolattice` (activate it, or prepend
its `bin/` to `PATH`); on another, system `python3` (with numpy) and `pdflatex` are available but conda is not.
Check with `which cmake python3 pdflatex` first.

All commands start from the repository root:
```bash
mkdir -p build_attractorE && cd build_attractorE && cmake .. -DMODEL=attractorE && make -j8 && cd ..
cd runs/E_starobinsky_gw
OMP_NUM_THREADS=8 OMP_PROC_BIND=false ../../build_attractorE/attractorE input=run.in overwriteFiles=true
cd ..   # now in runs/
python3 make_notebook.py && jupyter nbconvert --to notebook --execute --inplace attractor_runs.ipynb
cd ../notes && pdflatex implementation_report.tex && pdflatex implementation_report.tex   # rebuild a note (run twice for refs/TOC)
```
`make_notebook.py` writes the notebook next to itself, and the notebook reads the run folders relative to `runs/`, so both work from any clone location.
- **Performance** (first computer: 16 logical / 8 physical cores): **threads = physical cores is fastest**; using all logical cores is slower. **Never run two simulations at once**: they fight for cores and slow down by about 10×. N=64 with GWs takes ≈1.7 s per program-time unit, so ≈8.5 min to t=300.
- CosmoLattice refuses to overwrite output unless `overwriteFiles=true`.
- Any `.in` key can be overridden on the command line (`N=128 tMax=500`).
- Output: `average_energies.txt`, `average_energies_gws.txt` (rhoGW_over_rho, rhoGW), `average_scale_factor.txt`, and `spectra_scalar_{0,1}.txt` / `spectra_energy_gws.txt`. Spectra are blocks separated by blank lines with columns k, Δ_field, Δ_momentum, multiplicity; the GW spectrum columns are k, Ω_GW, multiplicity. Spectrum times are in `average_spectra_times.txt`.

## Shell gotchas found this session
- `pkill -f <pattern>` also kills the agent's own shell when the pattern appears in the command line (exit 144). Kill by PID: `ps -eo pid,args | awk '$2 ~ /build_attractor[ET]\/attractor[ET]$/{print $1}' | xargs -r kill`.
- `rm` on `$VAR/*` is blocked by a safety check; use literal paths or `"${VAR:?}"/*`.
- Foreground `sleep` is blocked; use background commands or Monitor.

## Key findings (do not rediscover)
1. **κ typo in 2510.18656:** the text says κ=0.9999 for all deformed benchmarks. For k=4,6 that gives n_s=0.99–1.04. **κ=0.99999 reproduces the quoted n_s, r**, and the figure file is named `fulltplots_kappa0p99999.pdf`. Exception: deformed E k=6, α=4 matches κ=0.9999. The paper's r≈0.012 for deformed T k=6 is not reproduced (we get 0.004).
2. **2π in 1707.04533 Eq. (2.17):** the numerical shortcut `f ≃ ε_i^{1/4}(k/ρ_i^{1/4})×8·10⁹ Hz` is 2π below the equation's own first line. The notebook uses the exact first line, which reproduces the standard Dufaux et al. ≈4·10¹⁰ Hz.
3. The helper reproduces the paper's Starobinsky numbers (λ=1.57e-10 vs 1.6e-10, x_*=5.32 vs 5.35 for N_*=55). Its n_s, r are slow-roll, accurate to ~1e-3.
4. Run results (N=64, k_IR=0.5, q₀≈1e4: E g=4.5e-3, T g=1.7e-3):
   - Broad resonance until t̃≈60–70, then backreaction; w reaches 0.25–0.27.
   - Final ρ_GW/ρ: 7.8e-5 (E) and 1.0e-4 (T).
   - Today, with instant RD after the run: peak f ≈ 2.6e9 Hz with h²Ω ≈ 1.6e-9 (E); peak f ≈ 1.35e9 Hz with h²Ω ≈ 2.1e-9 (T).

## Known limitations / open tasks
- **Runs are UV-limited:** spectra pile up near k_max ≈ 28 ω*, and the GW peak (k≈17) is close to it. Rerun at **N=128–256 with the same k_IR** before quoting numbers (N=128 ≈ 8× slower).
- Today's GW values assume instant radiation domination after t_f=290 (upper bound, since w≈0.24 then). The notebook has `N_POST`, `W_POST` knobs.
- N_* in the benchmark files was taken from the paper's reheating ranges. For k≥6 it should be fixed from the fragmentation temperature measured on the lattice.
- Benchmark batch (2026-10-03, see notebook Sec. 6): deformed k=2 and all k=4 points run with g (q0=1e4). k=4 needs dt=0.005 (χ frequency ≈2√q=200 in program time all run long; dt=0.01 is at the leapfrog limit). Friedmann violation ≲9e-3 for k=4 (∝dt²).
- **k≥6 at the paper's α=1–5 do not fragment**: g=0 self-resonance gives nothing to a≈34; with g, the real oscillation period is ≈20 program units (not 2π), q_eff is small and only k̃≲0.5 is unstable; a mode scan gives n growth ≤1e5 by t=150 even for g×100 (would need dt~1e-4). w stays at (k−2)/(k+2). Options: small α (strong self-resonance) or long g=0 baselines.
- The python with scipy/jupyter on the second computer is `~/anaconda3/envs/cosmolattice/bin/python` (system python3 has numpy only).
- Not yet run: `g` scans, and validation of `attractorT` (k=2, κ=1) against the existing `models/tanh2.h` (M=√(6α)MP, Λ⁴=1.5λMP⁴).
- A PR to upstream should probably exclude `runs/` and `papers/` (data and third-party sources).
- The user communicates in English, is a physicist (UBA), and prefers that the remote takes priority on pulls.
