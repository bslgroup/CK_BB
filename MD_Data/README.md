# Creatine Kinase (CK) — Molecular Dynamics of the WT Dimer

This repository contains the curated inputs, analysis scripts, and analysis
outputs for all-atom molecular dynamics (MD) simulations of the **WT creatine
kinase (CK-BB) homodimer**. It is a self-contained, reproducible package — everything needed to regenerate the
reported analyses **except the raw trajectories**, which are omitted for size
(see [Trajectories](#trajectories)).

## Systems

Twelve production runs: a 2 × 2 design (ligand state × pH), three independent
replicates each.

| System     | Ligand state                         | pH  | Protonation            | Replicates |
|------------|--------------------------------------|-----|------------------------|:----------:|
| `apo_5.5`  | apo (no substrate)                   | 5.5 | PROPKA set (see below) | 3          |
| `apo_7.3`  | apo (no substrate)                   | 7.3 | standard               | 3          |
| `crn_5.5`  | holo — ATP + creatine (CRN) + Mg²⁺   | 5.5 | PROPKA set (see below) | 3          |
| `crn_7.3`  | holo — ATP + creatine (CRN) + Mg²⁺   | 7.3 | standard               | 3          |

- **apo** = no bound substrate; **crn** = substrate-bound (ATP + creatine + Mg²⁺).
- **Protonation** was assigned with **PROPKA**, symmetric on both protomers. At
  **pH 5.5**, six titratable sites per chain are protonated —
  **Asp62, Glu19, Glu80, Glu155, Glu231, and His66** (structural numbering; add 3
  for the PSF numbering, which carries a 3-residue N-terminal `GLY-SER-HIS` tag
  remnant, so the native chain begins at Met = PSF resid 4). At **pH 7.3** none
  are protonated (standard states). Note: the non-native N-terminal tag histidine
  (PSF resid 3) is also protonated at pH 5.5 in the built system.

## Simulation details

- **Force field:** CHARMM36m; system built with CHARMM-GUI (`step3_input`), 310 K.
- **Engine:** **equilibration and the initial production were run with NAMD** on
  the Frontera supercomputer (TACC; job script `run_fr`); the remainder of the
  production was run on the **Anton 3** special-purpose supercomputer
  (D. E. Shaw Research; config `base.ark`).
- **Run inputs shipped** (`input/`): the NAMD `step4_equilibration.inp` and
  `step5.1.production.inp`, the full CHARMM36m `toppar/` set (crn systems also
  carry the ATP/creatine ligand parameters — `toppar/crn.prm` + `toppar/atp/`),
  the Frontera NAMD job script `run_fr`, and the Anton 3 config `base.ark`.
- **Structures shipped** (`input/`): the solvated build (`step3_input.psf/.pdb`)
  and the water-stripped pair used for analysis (`step3_input_nowat.psf/.pdb`).
- **Analysis trajectories:** PBC-wrapped and sub-sampled to **421 frames per
  replicate** (`dcd/ck_wrapped.dcd` in the working tree; not shipped here).

## Repository layout

```
Github/
├── README.md
└── <system>/                 # apo_5.5, apo_7.3, crn_5.5, crn_7.3
    ├── analysis/             # pipeline scripts (run from inside a replicate dir)
    ├── data/
    │   ├── 1/  2/  3/        # curated analysis outputs, one dir per replicate
    │   │   ├── rmsd.txt
    │   │   ├── rmsf_a.txt  rmsf_b.txt
    │   │   ├── EC-angle.txt
    │   │   ├── bridges/      # salt-bridge distances (.dat)
    │   │   └── hbond/        # crn systems only: ATP-phosphate↔His distances
    └── input/                # structures + simulation run inputs
        ├── step3_input(.psf/.pdb)         # solvated build
        ├── step3_input_nowat(.psf/.pdb)   # water-stripped (analysis)
        ├── step4_equilibration.inp        # NAMD equilibration
        ├── step5.1.production.inp         # NAMD production
        ├── toppar/                        # CHARMM36m params (+crn.prm & atp/ for crn)
        ├── run_fr                         # NAMD job script (Frontera/TACC)
        └── base.ark                       # Anton 3 run config
```

## Analyses

All analyses are VMD/Tcl driven by shell wrappers, with multi-replicate plotting
in gnuplot. Scripts live in each system's `analysis/`.

| Analysis        | Scripts (`analysis/`)                                   | Output (`data/<rep>/`)              |
|-----------------|---------------------------------------------------------|-------------------------------------|
| RMSD            | `rmsd.sh`, `RMSD.tcl`                                    | `rmsd.txt`                          |
| RMSF (per chain)| `rmsf.sh`, `RMSF_Apro.tcl`, `RMSF_Bpro.tcl`             | `rmsf_a.txt`, `rmsf_b.txt`          |
| Domain angle    | `angle.sh`, `EC-angle.tcl`, `macros_angle.tcl`          | `EC-angle.txt`                      |
| Salt bridges    | `salt_all.sh`                                            | `bridges/salt_*.dat`                |
| ATP–His contacts (crn only) | `hbond.sh`                                  | `hbond/dist_*.dat`                  |

- **RMSD** — whole dimer, Cα, per-chain (PROA/PROB) Cα, and per-chain helices,
  after alignment to frame 0.
- **RMSF** — per-residue Cα fluctuations, reported separately for the two
  protomers (chains PROA and PROB).
- **EC-angle** — inter-domain angle computed from the dot product of the inertia
  tensors of two "EC" regions.
- **Salt bridges** — minimum donor–acceptor distances for the
  **Asp329 – Ser202 – His69** network (files `salt_329-202`, `salt_202-69`,
  `salt_329-69`).
- **ATP-phosphate ↔ His194** (crn only) — distances from the ATP γ-phosphate
  oxygen (`O2p`) to His194 (`ND1`/`NE2`) on each protomer.

Shared VMD selection macros are in `macros_ck.tcl`; `remove_water.tcl`,
`z-wrap-pbc.tcl`, and `z-recenter.tcl` are the trajectory pre-processing
(water-strip / PBC-wrap / recenter) steps.

## Reproducing the analyses

The raw trajectories are not distributed here. To reproduce from a trajectory,
place `step3_input_nowat.psf` and a wrapped trajectory at `dcd/ck_wrapped.dcd`
inside a replicate directory alongside the `analysis/` scripts, then run e.g.:

```bash
bash rmsd.sh        # -> rmsd.txt
bash rmsf.sh        # -> rmsf_a.txt, rmsf_b.txt
bash angle.sh       # -> EC-angle.txt
bash salt_all.sh    # -> bridges/salt_*.dat  (+ multi-replicate plots)
bash hbond.sh       # crn systems: hbond/dist_*.dat
```

The `gnuplot-multi-{rmsd,rmsf,angle}.sh` wrappers regenerate the
multi-replicate comparison figures from the three `data/<rep>/` output sets.
