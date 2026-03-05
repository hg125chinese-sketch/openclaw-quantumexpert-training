# Aspirin (acetylsalicylic acid) — DFT end-to-end run (QE)

Molecule:
- name: aspirin (acetylsalicylic acid)
- SMILES: `CC(=O)Oc1ccccc1C(=O)O`

Method label (as executed here):
- DFT functional: **PBE0** (this run)  
  - Note: skill prefers D3(BJ) (or an XC with built-in dispersion). This environment lacked dftd3 at run time, so this archived run is **no-dispersion** and should be treated as a tooling-limited baseline.
- Basis sets: **def2-SVP** (DZ) and **def2-TZVP** (TZ)
- Dispersion: **none** (D3/BJ unavailable in this environment at run time)
  - Recommended rerun: **B97-D** or **WB97M-V** (built-in dispersion in XC) if you cannot enable D3.
- Software: **PySCF 2.12.1**
- SCF settings: grid level=4, conv_tol=1e-9

## Files
- `aspirin_dft_pipeline.py`
  - SMILES → RDKit 3D → PySCF+geomeTRIC geometry optimizations (SVP then TZVP)
  - property extraction (E, HOMO/LUMO/gap, dipole)
  - (initial version wrote summary, but ZPE required separate run; see below)

- `aspirin_freq_zpe.py`
  - recomputes harmonic frequencies + ZPE at the optimized geometries
  - updates `summary.json` to include `freqs_cm1`, `zpe_ha`, `zpe_kcal_mol`

- `aspirin_rdkit_init.xyz` — initial 3D structure from RDKit
- `aspirin_prod_def2-svp_opt.xyz` — optimized geometry at PBE0/def2-SVP
- `aspirin_prod_def2-tzvp_opt.xyz` — optimized geometry at PBE0/def2-TZVP
- `summary.json` — consolidated results for both basis sets

## Re-run (copy/paste)
From workspace root:

```bash
cd <repo_root>
python research/aspirin_dft/aspirin_dft_pipeline.py
python research/aspirin_dft/aspirin_freq_zpe.py
```

Notes:
- The minimum-frequency modes are low (~tens of cm^-1), so RRHO thermochemistry may be sensitive; this run reports ZPE only.
- Absolute energies between different basis sets are not directly comparable as a “convergence metric”; use relative energies on a fixed geometry if needed.
