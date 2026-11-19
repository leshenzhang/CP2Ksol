<p align="center"><img src="docs/logo/logo.png" alt="CP2Ksol" width="360"></p>

# CP2Ksol

[English](README.md) · [中文](README.zh.md)

Developed by [Leshen Zhang](https://github.com/leshenzhang) (Alexandrova group, UCLA).

Implicit solvent, implicit electrolyte and a constant-potential (grand-canonical) electrode for **CP2K 2026.2** (Quickstep, GPW). It runs single points, geometry optimisations and molecular dynamics, with or without explicit water in the cell.

The code is plain Fortran in the CP2K source tree. It needs no external library or embedded Python and leaves the MPI lifecycle unchanged.

Contents: [Features](#features) · [Install](#install) · [Quick start](#quick-start) · [Choosing a model](#choosing-a-model) · [Keywords](#keywords) · [Recommended settings](#recommended-settings) · [Performance and accuracy](#performance-and-accuracy) · [Output](#output) · [Tests](#tests) · [Known limitations](#known-limitations) · [Citing and license](#citing-and-license)

## Features

- **Implicit solvent.** Two cavity families. `METHOD ISODENSITY` puts one cavity at an isovalue of the electron density (https://doi.org/10.1063/1.4865107). `METHOD SOLVENT_EXCLUDED` builds layered closest-approach surfaces (solvent, dielectric, ion, cavitation) by probe-shell convolution (https://doi.org/10.1063/5.0176308).
- **Implicit electrolyte.** Linearised Poisson–Boltzmann screening with a G = 0 potential reference, so the bulk electrolyte potential is zero (https://doi.org/10.1063/1.5132354); size-modified nonlinear Poisson–Boltzmann ions (https://doi.org/10.1103/physrevlett.79.435); a nonlinear dielectric with dipolar saturation (https://doi.org/10.1088/0965-0393/21/7/074005).
- **Explicit water plus implicit solvent.** `&REGION` confines the implicit solvent to a slab; `&REGION &FLEXIBLE` makes the slab boundary follow the explicit water film; `&EXCLUSION` switches the implicit dielectric and ions off around explicit atoms that move. All three give analytic forces.
- **Constant potential.** `&DFT &SCF &GCE` makes the electron number fractional and steers the work function to a target potential. Schemes: `CAPACITANCE` (secant update of the electron number), `DYNAMICS` (electron number as a dynamical variable with a fictitious mass and an optional Nosé–Hoover chain). The electrolyte is the counter charge. The electron-number state is stored in restart files.
- **Solver.** `SOLVER NEWTON_CG`: Newton iterations on the true residual of the generalised Poisson equation, inner preconditioned CG, line search on the free energy. Speed-ups: ψ form, two-level preconditioner, potential history, tolerance coupled to the SCF.
- **Pseudo-core cavity density** for GTH pseudopotentials (`ADD_CORE_DENSITY`), so the cavity cannot open inside an atom.

## Install

The code is 6 new modules and 23 modified CP2K files, all under `src/` (`src/CMakeLists.txt` lists the new files). Start from an unpacked CP2K 2026.2 source tree.

```bash
tar xjf cp2k-2026.2.tar.bz2
tools/install.sh cp2k-2026.2          # copies src/ over the matching files, keeps *.orig backups
# or, equivalently:
cd cp2k-2026.2 && patch -p1 < ../patches/cp2k-2026.2.diff
```

`install.sh` refuses any tree whose `src/cp2k_info.F` is not version 2026.2. Then build CP2K as usual (CMake, see `INSTALL.md` of CP2K) and use the `psmp` build (MPI and OpenMP).

Toolchains the code has been built and run with:

| Compiler and MPI | Libraries | Flags |
|---|---|---|
| gcc 11.4, Intel MPI 2021.7.1 | MKL 2022.2.1, FFTW 3.3.10 | `-O3 -march=x86-64-v3 -mtune=generic` |
| gcc 13, OpenMPI 4.1 | OpenBLAS, build type GENERIC | `-O3 -march=x86-64-v3 -mtune=generic` |
| Cray MPICH | Cray libraries | `-O3 -march=znver3` |

On Cray MPICH the program sets `MPICH_SMP_SINGLE_COPY_MODE=NONE` before MPI starts when the variable is unset, `CMA` or `XPMEM`. The single-copy transfers leak memory in ScaLAPACK diagonalisations; after about 4 GiB are lost the ranks wait in BLACS for ever, which stalls long constant-potential runs. Any other value set in the job script is kept.

Layouts are job settings; the code defaults do not depend on them. Results of other layouts agree with the 64 × 1 default within |dE| ≤ 3.5e-10 Ha and max |dF| ≤ 8.8e-7 Ha/bohr. For a 64-core node, 32 ranks × 2 threads with `OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close` and `MAX_DISTRIBUTED_LEVEL 1` under `&DFT &MGRID &RS_GRID` is a good start; a 16-core node runs 16 × 1 at the defaults.

To regenerate the patch after editing `src/`, run `tools/make_patch.sh <pristine cp2k-2026.2 root>`.

## Quick start

Examples 1 and 2 are the `&SCCS` (and `&SCF`) part of a deck in `tests/`; cell, coordinates, basis sets and the rest of the input stay as in that file. Example 3 is the production setting of a 203-atom film. All three parse with `cp2k.psmp --check`.

### 1. A molecule in implicit solvent

Base: `tests/regtest-sccs-lpb-1/H2O_isodensity_lpb.inp` (H2O in a 10 Å box, PBE, `CUTOFF 400`). Replace its `&SCCS` block by:

```
&SCCS on
  METHOD ISODENSITY
  RELATIVE_PERMITTIVITY 78.4
  GAMMA [mN/m] 18.47
&END SCCS
```

This is the ISO model. `SOLVER NEWTON_CG` is the default for `ISODENSITY`; `GAMMA` is the cavitation surface tension fitted for solvation free energies of molecules. The run converges in 18 SCF steps to about −17.2305 Ha.

### 2. A charged metal surface at constant potential

Base: `tests/regtest-sccs-lpb-3/Al111_gce_cap_lpb.inp` (Al(111) 2×2, 3 layers, Γ point). Replace its `&SCCS` and `&SCF` blocks by the following. A constant-potential run needs `ADDED_MOS`, `&DIAGONALIZATION` and `&SMEAR METHOD FERMI_DIRAC`:

```
&SCCS on
  METHOD ISODENSITY
  RELATIVE_PERMITTIVITY 78.4
  DEBYE_LENGTH [angstrom] 3.0       ! 1 M 1:1 electrolyte, the counter charge of the electrode
&END SCCS
&SCF
  EPS_SCF 1.0E-6
  MAX_SCF 150
  ADDED_MOS 20
  &SMEAR on
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 300
  &END SMEAR
  &DIAGONALIZATION on
    ALGORITHM STANDARD
  &END DIAGONALIZATION
  &MIXING
    METHOD BROYDEN_MIXING
    ALPHA 0.2
    NBROYDEN 8
  &END MIXING
  &GCE on
    UPDATE_SCHEME CAPACITANCE
    UPDATE_EVERY SCF_CYCLE
    UPDATE_INTERVAL 5
    TARGET_POTENTIAL [eV] -0.6      ! U = -0.6 V vs SHE
    PHI_SHE [eV] 4.600
    FERMI_TOL [eV] 0.01
    MAX_OUTER_UPDATES 3
  &END GCE
&END SCF
```

The electron number N_e is adjusted until the work function matches the target. Each ionic step prints a `GCE| trace` line with the potential reached ([Output](#output)); the run ends with |U − U_T| below 0.01 V.

### 3. An explicit water film at constant potential (SE with explicit water)

This is the setting used for the 203-atom Cu(111)/water/K⁺ film in the paper: a 4×4×4 Cu slab (the lowest 32 atoms fixed), 46 water molecules and one K⁺ in a 10.267 × 8.892 × 47.000 Å cell, Γ point, PBE with DFT-D3(BJ), `CUTOFF 400`, `REL_CUTOFF 40`, NVT at 300 K with a CSVR thermostat, 1 fs steps, held at −0.600 V vs SHE with the electron number updated every 5 steps. Kinds are `Cu`, `K`, `Ow` (water O) and `Hw` (water H). Only the implicit-solvent, `&EXCLUSION`, `&REGION`/`&FLEXIBLE` and constant-potential blocks are given; cell, coordinates, `&MOTION` and the kind definitions are those of an ordinary MD input.

```
&SCCS on
  RELATIVE_PERMITTIVITY 78.4
  GAMMA [mN/m] 0.0
  DELTA_RHO 2.0E-5
  DERIVATIVE_METHOD fft
  EPS_SCCS 1.0E-6
  EPS_SCF 0.5
  MAX_ITER 200
  SOLVER NEWTON_CG
  SOLVER_TOLERANCE 1.0E-8
  SOLVER_MAX_ITER 50
  DEBYE_LENGTH [angstrom] 3.0
  NONLINEAR_IONS F
  NONLINEAR_DIELECTRIC F
  SOLVER_GRID_LEVEL 1
  METHOD SOLVENT_EXCLUDED
  &SOLVENT_EXCLUDED
    RHO_CAVITY [angstrom^-3] 0.0158
    SWITCH_WIDTH 0.6
    PROBE_KERNEL_WIDTH [angstrom] 0.125
    SOLVENT_RADIUS [angstrom] 1.4
    ION_RADIUS [angstrom] 4.0
    DIELECTRIC_RADIUS [angstrom] 0.7
    CAVITATION_RADIUS [angstrom] 0.0
  &END SOLVENT_EXCLUDED
  &REGION
    AXIS 3
    WIDTH [angstrom] 0.5
    PROBE_ATOMS 65                  ! the K+ ion; reporting only
    PROBE_RADIUS [angstrom] 3.0
    &FLEXIBLE
      KIND Ow 1.0                   ! solvent oxygens and ions only
      KIND K 1.0
    &END FLEXIBLE
  &END REGION
  &EXCLUSION on
    KIND Ow 1.80 2.40               ! kind, r_diel, r_ion in angstrom
    KIND Hw 1.10 1.70
    KIND K  2.00 2.60
    WIDTH [angstrom] 0.25
    PRINT_SUMMARY T
  &END EXCLUSION
  SCF_TOL_FACTOR 0.1
  SOLVER_STRIDE 1
  PSI_FORM T
&END SCCS
&SCF
  &GCE on
    UPDATE_SCHEME CAPACITANCE
    UPDATE_EVERY SCF_CYCLE
    UPDATE_INTERVAL 5               ! update N_e every 5 MD steps
    TARGET_POTENTIAL [eV] -0.600    ! U = -0.6 V vs SHE
    PHI_SHE [eV] 4.600
    CAPACITANCE_INIT 0.4749
    FERMI_TOL [eV] 0.01
    FERMI_UPDATE_TOL [eV] 0.005
    MIN_STEPS_AFTER_UPDATE 12
    MAX_OUTER_UPDATES 6
    EXCESS_ELECTRONS 0.673537       ! starting N_e of this frame
  &END GCE
  EPS_SCF 1.0E-6
  MAX_SCF 200
  SCF_GUESS ATOMIC
  ADDED_MOS 117
  &SMEAR on
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 300
  &END SMEAR
  &DIAGONALIZATION
    ALGORITHM STANDARD
  &END DIAGONALIZATION
  &MIXING on
    METHOD BROYDEN_MIXING
    ALPHA 0.2
    BETA 1.5
    NBROYDEN 8
  &END MIXING
&END SCF
```

The solvent is SE with `DIELECTRIC_RADIUS` 0.7, and `&EXCLUSION` and `&REGION &FLEXIBLE` add the exclusion spheres and the atom-following boundary for the explicit water. `&FLEXIBLE` lists the water oxygens and the K⁺ ion and leaves out the Cu atoms. `EXCESS_ELECTRONS` and `CAPACITANCE_INIT` are starting values for this frame; a new system needs its own, or the defaults. The complete input, with cell and coordinates, is in `examples/film203_se/`. For a smaller test of the same blocks use `tests/regtest-sccs-lpb-5/K6H2O_dbg_se_exclusion_linear_all.inp` (K⁺ and 6 H2O) as the base.

## Choosing a model

| Model | Keywords | Parameters |
|---|---|---|
| ISO | `METHOD ISODENSITY`, linear dielectric, linear Poisson–Boltzmann | n_c = 3.00e-3 Å⁻³ (`RHO_CAVITY`, default), σ = 0.6 (`SWITCH_WIDTH`) |
| SE | `METHOD SOLVENT_EXCLUDED`, linear dielectric (`NONLINEAR_DIELECTRIC F`), `DIELECTRIC_RADIUS 0.7`, linear Poisson–Boltzmann | n_c = 1.58e-2 Å⁻³, σ = 0.6, R_solv = 1.4 Å (`SOLVENT_RADIUS`), R_ion = 4.0 Å (`ION_RADIUS`), R_diel = 0.70 Å. With explicit solvent atoms, add `&EXCLUSION` and `&REGION &FLEXIBLE` (exclusion spheres and atom-following boundary); `&FLEXIBLE` lists only the solvent oxygens and the ions. Without explicit solvent atoms this is plain SE |
| SE-NLD | `METHOD SOLVENT_EXCLUDED`, nonlinear dielectric (the default of this method) | R_diel = 1.0 Å (default) |
| SE-NLPB | SE-NLD with `NONLINEAR_IONS T`, `ION_STERIC_RADIUS 4.0` | size-modified nonlinear ions |

The bulk permittivity is ε_b = 78.4 for all models (`RELATIVE_PERMITTIVITY 78.4`).

One rule covers most cases: ISO when the cell has no explicit solvent, SE (with `&EXCLUSION` and `&REGION &FLEXIBLE`) when it has explicit water.

- Interfaces and molecules without explicit solvent: ISO.
- Cells with explicit water: SE with `&EXCLUSION` and `&REGION &FLEXIBLE`. List only the solvent oxygens and the ions in `&FLEXIBLE`, not the electrode atoms.
- Isolated cations: ISO with `RHO_CAVITY [angstrom^-3] 0.0010`. SE is too strong for isolated small cations.
- SE-NLD and SE-NLPB add dipolar saturation and finite-size ions. On Ag(100) the linear models give a flat C(U), SE-NLD a bell-shaped one, and SE-NLPB a dip at the potential of zero charge.

## Keywords

The full reference is in [`docs/keywords.md`](docs/keywords.md) ([中文](docs/keywords.zh.md)): every keyword with meaning, default and unit, the former keyword names, and a symbol table. Sections:

| Section | What it controls |
|---|---|
| `&DFT &SCCS` | model, solver, electrolyte, nonlinear dielectric, pseudo-core cavity, speed-ups, diagnostics |
| `&SCCS &ISODENSITY`, `&SCCS &SOLVENT_EXCLUDED` | cavity parameters of the two methods |
| `&SCCS &REGION`, `&REGION &FLEXIBLE`, `&SCCS &EXCLUSION` | explicit water inside an implicit solvent |
| `&DFT &SCF &GCE` | constant potential: schemes, target, update loop, restart state, electron-number dynamics |

## Recommended settings

Everything not shown stays at its default. Water at 298 K: `RELATIVE_PERMITTIVITY 78.4`; 1 M 1:1 electrolyte: `DEBYE_LENGTH [angstrom] 3.0`. With an electrolyte use `&POISSON PERIODIC XYZ`, `POISSON_SOLVER PERIODIC`.

**Potential scale.** U = W − Φ_SHE, with W the work function and Φ_SHE = 4.600 eV, the default of `PHI_SHE`. `TARGET_POTENTIAL [eV]` takes U in V vs SHE.

**Surface tension.** Slabs, adsorption energies, potentials of zero charge and capacitances: `GAMMA` 0 (the default). Solvation free energies of molecules: `GAMMA [mN/m] 18.47` with ISO, `GAMMA [mN/m] 11.09` with SE.

**Constant-potential MD.** `UPDATE_SCHEME CAPACITANCE`, `UPDATE_EVERY SCF_CYCLE`, `UPDATE_INTERVAL 5`: the electron number is updated every 5 MD steps. The code default is 1 (update at every step); 5 is the recommended value, a balance between speed and accuracy. On a 203-atom Cu/water/K⁺ film at −0.600 V (2 ps, 1 fs steps) updating at every step holds the potential with a fluctuation of 5.10 mV at 42.18 s per step; updating every 5 steps gives 22.40 mV at 32.54 s per step, 22.9 % less time. `MAX_SCF` counts the inner steps of one update cycle, so 60 is enough. A restart file carries the whole state, so a restartable `GEO_OPT` or MD needs no extra keyword.

For single points on large cells that start far from the target, add to the `&GCE` block:

```
PIN_MIXER_G0 T
EARLY_UPDATE_DELTA 3.0E-3
EARLY_UPDATE_MIN_STEPS 6
EARLY_UPDATE_WF_STABLE [eV] 0.003
EARLY_UPDATE_MAX 4
EARLY_UPDATE_SECANT T
CAP_REFRESH_MIN [eV] 0.05
FINAL_TOL_FRAC 0.25
```

On a 323-atom Cu/water/K⁺ film (ISO, 42 cores) this cut the SCF steps from 164 to 92 and the wall time from 862 to 533 s. Cells that start within 0.02 e of the target gain 0–23 % of the steps and can be slower, so compare on one point first.

**Electron-number dynamics.** `UPDATE_SCHEME DYNAMICS` with `POTENTIOSTAT_PERIOD` (well above 10 time steps), `POTENTIOSTAT_CHAIN_LENGTH 3` for canonical sampling and `POTENTIOSTAT_CAPACITANCE` from a static N_e scan.

**Not recommended.** `UPDATE_SCHEME MU_MIXING` with the electrolyte as counter charge: the electron number moves in every SCF step, and a transient value can leave no volume for the ions, which stops the run. Use `CAPACITANCE`.

## Performance and accuracy

| Quantity | Result |
|---|---|
| Solvation free energy, FreeSolv, 222 molecules, ISO | MAE 1.52 kcal/mol |
| Potential of zero charge, nine Cu, Ag and Au faces against experiment, ISO, Φ_SHE = 4.600 eV | RMS 0.090 V |
| Continuum solve, Pt(111) 3×3, 3 layers, ISO with 1 M electrolyte | 0.128 s per SCF step |
| Constant-potential MD, 203-atom film, update every step / every 5 steps | 42.18 s per step, σ(U) 5.10 mV / 32.54 s per step, σ(U) 22.40 mV |

## Output

One line per ionic step records the potential:

```
GCE| trace  step <n>  N_e <x> e  W_f <y> eV  U <z> V vs SHE  U_target <t>  dU <d> V
GCE| omega  step <n>  E <e> Ha  n_exc <x> e  Omega <w> Ha  WF_target <t> eV
GCE| conserved  MD step <m>  Kin <k> Ha  E <e> Ha  dE <d> Ha  n_exc <x> e  Omega <w> Ha  Cons <c> Ha  Cons_Omega <co> Ha
GCE| dyn  step <n>  md_step <i>  N_e <x> e  dU <y> V  meas <T|F>  P_n <p>  Omega <o> Ha  Hextra <h> Ha  K_n <k>  E_th <t>
```

`trace` and `omega` are printed once per ionic step by every run with `&GCE` on, `conserved` by MD runs, `dyn` by `DYNAMICS`. `GCE| outer update`, `GCE| early update` and `GCE| carry` mark N_e updates; `GCE| state read from the restart` is printed after a restart.

`Omega = E + WF_target n_exc` with `n_exc = N_e − N_input` and `WF_target = −mu_e`. With a varying N_e the conserved quantity is `Cons_Omega = Cons + WF_target n_exc`, not `Cons`.

## Tests

`tests/regtest-sccs-lpb-{1..7}` follow the CP2K regression-test layout and are read through `TEST_FILES.toml`: an energy tier on H2O, finite-difference force tiers (`RUN_TYPE DEBUG`), constant-potential tiers on Al(111), guards, the exclusion mask, the flexible region and the restart state. Copy them to `tests/QS/` of the CP2K tree, add `QS/regtest-sccs-lpb-1` … `QS/regtest-sccs-lpb-7` to `tests/TEST_DIRS`, and run the CP2K regression tester as usual.

Decks that must abort, and restart decks that need the restart file of an earlier deck, are not listed in `TEST_FILES.toml`; run them by hand as described in `tests/README.md`, which also lists the tolerances and the conditions that a matcher cannot express (for example `grep -c "GCE| outer update"` must be > 0 for the Al(111) deck).

`tests/fivepoint/` holds the five-point charging scan of Pt(111) 4×4×4 with SE-NLPB: `cp2k.psmp -i fp_x00.inp -o fp_x00.out` for each deck, then `python3 fit_fivepoint.py fp_*.out`. It reproduces C = 1.1764 e/V = 17.70 µF/cm² and U_pzc = 0.8500 V vs SHE.


## Citing and license

Methods the implementation builds on:

- Isodensity cavity, pseudo-core density: https://doi.org/10.1063/1.4865107
- Implicit electrolyte, linearised Poisson–Boltzmann with a G = 0 potential reference: https://doi.org/10.1063/1.5132354
- Solvent-excluded cavities, nonlinear dielectric, size-modified ions: https://doi.org/10.1063/5.0176308
- Nonlinear fluid response: https://doi.org/10.1088/0965-0393/21/7/074005
- Size-modified Poisson–Boltzmann equation: https://doi.org/10.1103/physrevlett.79.435
- Solution scheme for the size-modified Poisson–Boltzmann equation: https://doi.org/10.1021/acs.jctc.6b00435
- Self-consistent continuum solvation, `METHOD ANDREUSSI`: https://doi.org/10.1063/1.3676407
- Generalised Poisson and Poisson–Boltzmann solver, ψ form: https://doi.org/10.1063/1.4939125
- Grand-canonical ensemble approaches in CP2K: https://doi.org/10.1021/acs.jctc.4c00671
- Solvent-aware SCCS in CP2K, `METHOD SAA_ANDREUSSI`: https://doi.org/10.1016/j.cpc.2025.109563
- Nosé–Hoover chain integrator: https://doi.org/10.1080/00268979600100761

CP2Ksol is developed and maintained by Leshen Zhang (https://github.com/leshenzhang). If you use this code, cite this repository (see `CITATION.cff`) and the papers above that apply to the models you run.

GPL-3.0-or-later. CP2Ksol modifies CP2K, which the CP2K developers distribute under GPL-2.0-or-later, a license that permits redistribution under version 3 or any later version. The modified CP2K files keep their upstream headers, and the modifications and the new modules in this repository are licensed under GPL-3.0-or-later. See `LICENSE`.
