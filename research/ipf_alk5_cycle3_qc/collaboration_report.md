# IPF/ALK5 Cycle 3 — QE↔CE Collaboration Report (Batch QC Screening)

**Date (UTC):** 2026-03-04  
**Project:** `IPF_ALK5_cycle3`  
**Upstream:** ChemicalExpert (CE) VAE-generated candidates (RDKit-preprocessed, still potentially problematic)  

## 0) Deliverables (what was produced)

Written to `exports/` (handoff artifacts for CE):
- `exports/qc_results.json` — per-molecule QC results (QE→CE data contract)
- `exports/batch_audit.json` — audit log with per-step events and failure reasons

Per-molecule checkpoints (for restart / forensic debugging):
- `exports/checkpoints_cycle3/<mol_id>.json`
- `exports/checkpoints_cycle3/<mol_id>_init.xyz`
- (when optimization succeeded) `exports/checkpoints_cycle3/<mol_id>_opt.xyz`

Source input file (from CE):
- `imports/qc_handoff_cycle3.json`

## 1) Workflow executed (qchem-workflow: batch screening pipeline)

We executed a checkpointed batch pipeline consistent with `qchem-workflow` “Batch QC screening” philosophy:

1) **Input validation / 3D generation**: RDKit parse + ETKDG embedding (+ fail-soft MMFF/UFF cleanup)
2) **DFT optimization**: **B97-D/def2-SVP** (PySCF + geomeTRIC via PySCF wrapper)
3) **Frequency analysis**: harmonic Hessian at optimized geometry (gate: `n_imag == 0` for minima)
4) **Property extraction** (on optimized geometry): `E_total`, `HOMO`, `LUMO`, `gap`, `dipole`
5) **qc_flag per molecule** (batch never stops on a failure)
6) **Checkpoint strategy**: write after each molecule (results + audit), so partial progress is always saved

### Method label (applies to all PASS molecules)
- **XC:** B97-D (dispersion-inclusive XC)
- **Basis:** def2-SVP
- **Software:** PySCF 2.12.1, geomeTRIC 1.1

## 2) Batch summary

**N molecules received:** 4  
**N PASS (opt + freq gates):** 2  
**N OPT_FAIL:** 2  
**Observed OPT_FAIL rate:** 50%  

Failure modes seen (from audit):
- `cycle3_top5_hinge_1`: geomeTRIC/PySCF gradient gate failure due to unconverged gradient (`Nuclear gradients ... not converged`).
- `cycle3_top5_hinge_3`: SCF did not converge at initial geometry; optimization could not start.

## 3) Per-molecule results (QE→CE contract fields)

Values below are directly from `exports/qc_results.json`.

### PASS molecules

**(A) cycle3_top5_hinge_2**
- Input: charge=0, spin=0; vina_score=-8.003; hinge_hbond=true
- Method label: **B97-D/def2-SVP** (PySCF 2.12.1)
- Gates: opt=PASS; freq=PASS (`n_imag = 0`)
- Properties:
  - `E_total = -1277.2308967891386 Ha`
  - `HOMO = -0.16206881042624227 Ha`
  - `LUMO = -0.08529378601133561 Ha`
  - `gap = 2.089154843399184 eV`
  - `dipole = 1.3330793062611015 D`
  - `min_freq = 28.167441058142963 cm^-1`

**(B) cycle3_top5_hinge_4**
- Input: charge=0, spin=0; vina_score=-8.295; hinge_hbond=true
- Method label: **B97-D/def2-SVP** (PySCF 2.12.1)
- Gates: opt=PASS; freq=PASS (`n_imag = 0`)
- Properties:
  - `E_total = -1067.916170660858 Ha`
  - `HOMO = -0.16005195339340328 Ha`
  - `LUMO = -0.05835218022792183 Ha`
  - `gap = 2.7673918087352813 eV`
  - `dipole = 3.8046845982393647 D`
  - `min_freq = 11.340802216057801 cm^-1`

### OPT_FAIL molecules (skipped, batch continued)

**(C) cycle3_top5_hinge_1**
- qc_flag: `OPT_FAIL`
- Audit failure reason (batch_audit):
  - `RuntimeError('Nuclear gradients of <...RKS_Scanner...> not converged')`

**(D) cycle3_top5_hinge_3**
- qc_flag: `OPT_FAIL`
- Audit failure reason (batch_audit):
  - `RuntimeError('SCF did not converge at initial geometry; cannot start optimization')`

## 4) Collaboration test outcome (QE↔CE loop closure)

CE reported:
- QE’s QC results were **received and integrated** into the multi-objective scoring.
- `hinge_4` ranked **#1** after CE aggregation.
- Based on the observed **~50% OPT_FAIL rate**, CE added a **MACE-OFF geometry pre-screen gate** prior to the next handoff.

**End-to-end status:** The collaboration test achieved a full closed-loop run:  
CE generation → handoff JSON → QE QC batch (checkpointed) → CE ingestion → ranking update → pipeline refinement (pre-screen gate).

## 5) Notes / recommendations for next cycle

1) **Pre-screening is justified:** VAE molecules can have hard SCF/geometry pathologies even after RDKit cleanup. A MACE-OFF geometry gate should reduce DFT failures and wasted CPU.
2) For remaining OPT_FAIL cases, if CE wants salvage attempts, recommend a tiered retry policy (opt-in):
   - try alternative SCF stabilizers (level shift, damping, DIIS tweaks),
   - try different initial geometry (multiple ETKDG seeds),
   - fail-fast if SCF still not converged.
3) Keep writing per-molecule checkpoints; it materially reduces operational risk for long optimizations/frequencies.
