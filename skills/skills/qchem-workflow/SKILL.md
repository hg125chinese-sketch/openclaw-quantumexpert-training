---
name: qchem-workflow
description: Orchestrate multi-step quantum chemistry workflows (geometry → frequency → excited state → benchmarking), manage batch calculations with QC gates at every step, and define standardized data interfaces for cross-agent collaboration with ChemicalExpert.
homepage: https://pyscf.org
metadata: { "openclaw": { "emoji": "🔗", "requires": { "bins": ["python3"], "python": ["pyscf", "numpy", "pandas", "json"] } } }
---

# Quantum Chemistry Workflow Orchestration

The "conductor" skill — orchestrates multi-step quantum chemistry calculations, enforces gates between steps, manages batch processing, and defines the data contract between QuantumExpert and other agents (especially ChemicalExpert). This skill does not compute anything itself; it calls qchem-dft, qchem-excited-state, and vault skills in the right order with the right checks.

## When to Use

- User requests a multi-step QC pipeline (optimize → frequency → properties → excited states)
- User sends a batch of molecules for QC-level evaluation
- ChemicalExpert sends candidates for electronic structure screening
- User wants to compare multiple methods/basis sets systematically
- User needs a reproducible workflow with full audit trail
- Any task that chains 2+ QC calculations where output of step N feeds step N+1

## Core Philosophy

1. **Plan the full pipeline before computing anything.** List all steps, gates between them, expected outputs. Get user approval. Then execute.
2. **Gates between every step.** Never feed garbage forward. If step N fails its gate, step N+1 does not start.
3. **Batch = many single molecules, not one big calculation.** Each molecule is independent. One failure does not stop the batch — flag it and continue.
4. **Standardized output format.** Every result must be machine-readable (JSON) and human-readable (summary table). Both get written.
5. **Audit trail is mandatory.** Log: which skills called, in what order, with what parameters, what passed/failed, wall time. This is the lab notebook.
6. **Cross-agent data contracts are explicit.** When receiving from or sending to CE, the format is defined here — not improvised per task.

## Phase 1: Pipeline Templates

### 1.1 Standard Single-Molecule Pipeline

```
Standard QC evaluation pipeline:
│
├── Step 0: Input validation
│   ├── SMILES → RDKit parse → 3D embed (ETKDG)
│   ├── Gate: RDKit mol is not None; 3D embed succeeded
│   └── Output: initial_geometry.xyz, charge, spin
│
├── Step 1: Geometry optimization [qchem-dft]
│   ├── Method: functional/basis per qchem-dft decision tree
│   ├── Gate: geomeTRIC "Converged!"; SCF converged; reasonable bond lengths
│   └── Output: optimized_geometry.xyz
│
├── Step 2: Frequency analysis [qchem-dft]
│   ├── Method: same as Step 1
│   ├── Gate: n_imag == 0 (for minimum); SCF converged
│   ├── If n_imag > 0: follow qchem-dft Phase 5.2 (displace + re-optimize)
│   └── Output: frequencies, ZPE, thermochemistry
│
├── Step 3: Single-point properties [qchem-dft]
│   ├── Method: same or higher level (e.g., larger basis on Step 1 geometry)
│   ├── Gate: SCF converged; spin contamination < 0.1 (if open-shell)
│   └── Output: E_total, HOMO, LUMO, gap, dipole, partial charges
│
├── Step 4 (optional): Excited states [qchem-excited-state]
│   ├── Method: TD-DFT or CASSCF/NEVPT2 per decision tree
│   ├── Gate: all requested states converged; NTO analysis done
│   └── Output: excitation energies, oscillator strengths, state characters
│
├── Step 5 (optional): Basis set convergence [qchem-dft]
│   ├── Same geometry, SVP → TZVP → TZVPP single points
│   ├── Gate: Δ(TZ→QZ) within thresholds from qchem-dft Phase 2.3
│   └── Output: convergence table
│
└── Step 6: Summary + archive
    ├── Write: summary.json (machine-readable) + summary.md (human-readable)
    ├── Write: audit_log.json (skills called, parameters, pass/fail, timing)
    └── Git commit to research/ directory
```

### 1.2 Batch Screening Pipeline (for CE integration)

```
Batch QC screening (e.g., 50 molecules from CE):
│
├── Step 0: Receive input from CE
│   ├── Expected: CSV or JSON with columns:
│   │     - mol_id
│   │     - smiles
│   │     - charge
│   │     - spin
│   │     - (optional) geometry_source: "rdkit_embed" | "diffsbdd_3d"
│   │     - (optional) geometry_file: path to SDF/XYZ provided by CE
│   ├── Validate:
│   │     - all SMILES parseable (if provided)
│   │     - geometry_file exists and is readable (if provided)
│   │     - no duplicate mol_id
│   └── Output: validated_input.json (N molecules)
│
├── Step 1: Initial geometry (choose one)
│   ├── If CE provided 3D geometry (geometry_source + geometry_file):
│   │     - Gate: geometry can be loaded; reasonable bond lengths; no severe clashes
│   │     - Output: initial_geometry.xyz (from CE geometry)
│   │     - NOTE: skip RDKit embedding
│   └── Else (no 3D geometry provided):
│         - SMILES → RDKit embed (ETKDG)
│         - Gate: embed succeeded; no severe clashes
│         - Output: initial_geometry.xyz (RDKit)
│
├── Step 2: DFT optimization + frequency [qchem-dft]
│   ├── Method: B3LYP-D3(BJ)/def2-SVP (restored default now that pyscf-dispersion is verified in QE) or per project spec
│   ├── Batch strategy: process one molecule at a time
│   │   ├── If SCF fails: log, skip, continue
│   │   ├── If optimization doesn't converge: log, skip, continue
│   │   ├── If n_imag > 0: attempt one displacement fix; if still fails, log
│   │   └── Each success: write mol_id_opt.xyz + freq data
│   ├── Gate per molecule: SCF converged AND n_imag == 0
│   └── Output: optimized_geometries/, per-molecule results
│
├── Step 3: Property extraction [qchem-dft]
│   ├── For each optimized molecule: E_total, HOMO, LUMO, gap, dipole
│   ├── Gate: SCF converged at optimized geometry
│   └── Output: properties.csv
│
├── Step 4: QC summary + delivery to CE
│   ├── Merge all results into: qc_results.csv
│   │   Columns: mol_id, smiles, method, basis, software, version,
│   │            E_total_Ha, HOMO_Ha, LUMO_Ha, gap_eV, dipole_D,
│   │            ZPE_kcal, n_imag, converged, qc_flags
│   ├── qc_flags examples:
│   │   - "PASS" — all gates passed
│   │   - "SCF_FAIL" — SCF did not converge
│   │   - "IMAG_FREQ" — imaginary frequency found
│   │   - "SPIN_CONTAM" — spin contamination > threshold
│   │   - "OPT_FAIL" — geometry optimization did not converge
│   ├── Write: batch_audit.json (per-molecule timing, pass/fail, skip reason)
│   └── Send to CE: qc_results.csv + batch_audit.json
│
└── Step 5: Batch statistics
    ├── Report: N_total, N_passed, N_failed (by failure type)
    ├── Flag: if failure rate > 20%, something systemic may be wrong
    └── Output: batch_summary.md
```

### 1.3 Method Comparison Pipeline

```
Systematic method comparison (e.g., functional benchmark):
│
├── Step 0: Define comparison matrix
│   ├── Molecules: list of test molecules with experimental reference values
│   ├── Methods: list of functionals/basis sets to compare
│   ├── Property: what to compare (gap, barrier, geometry, excitation energy)
│   └── Output: comparison_plan.json
│
├── Step 1: Run all combinations
│   ├── For each (molecule, method) pair:
│   │   ├── Optimize geometry (if geometry-dependent property)
│   │   ├── Compute target property
│   │   ├── Gate: SCF converged, all QC checks pass
│   │   └── Log result + timing
│   └── Output: raw_results.csv
│
├── Step 2: Statistical analysis
│   ├── For each method:
│   │   ├── MAE (mean absolute error vs reference)
│   │   ├── RMSE
│   │   ├── MAX error
│   │   ├── Systematic bias (signed mean error)
│   │   └── Outliers (error > 2σ)
│   └── Output: benchmark_stats.csv
│
├── Step 3: Recommendation
│   ├── Best method for this property/system class
│   ├── Cost-accuracy tradeoff (wall time vs MAE)
│   └── Output: benchmark_report.md
│
└── Gate: all molecules × methods completed (or failures documented)
```

## Phase 2: Inter-Step Gates

### 2.1 Gate Definitions

Every gate is a binary pass/fail check between pipeline steps. If a gate fails, the molecule does not proceed to the next step.

```python
#!/opt/conda/envs/chem/bin/python
"""Standard gate checks for QC workflow pipelines."""

def gate_scf_converged(mf):
    """Gate: SCF must converge."""
    if not mf.converged:
        return False, "SCF_FAIL: not converged after max_cycle iterations"
    return True, "PASS"

def gate_geometry_converged(opt_log):
    """Gate: geometry optimization must converge."""
    # Check for geomeTRIC "Converged!" in log
    if "Converged!" not in opt_log:
        return False, "OPT_FAIL: geomeTRIC did not report convergence"
    return True, "PASS"

def gate_frequency(freqs, mode="minimum"):
    """Gate: frequency analysis.
    mode='minimum': n_imag must be 0
    mode='ts': n_imag must be exactly 1
    """
    n_imag = sum(1 for f in freqs if f < 0)
    if mode == "minimum" and n_imag != 0:
        return False, f"IMAG_FREQ: n_imag={n_imag} (expected 0 for minimum)"
    if mode == "ts" and n_imag != 1:
        return False, f"IMAG_FREQ: n_imag={n_imag} (expected 1 for TS)"
    return True, "PASS"

def gate_spin_contamination(mf, threshold=0.1):
    """Gate: spin contamination check for open-shell systems."""
    if mf.mol.spin == 0:
        return True, "PASS (closed-shell, no check needed)"
    s2, _ = mf.spin_square()
    s_ideal = mf.mol.spin / 2.0
    s2_ideal = s_ideal * (s_ideal + 1)
    deviation = abs(s2 - s2_ideal)
    if deviation > threshold:
        return False, f"SPIN_CONTAM: <S²>={s2:.4f}, expected {s2_ideal:.4f}, Δ={deviation:.4f}"
    return True, "PASS"

def gate_basis_convergence(energies_by_basis, threshold_kcal=0.5):
    """Gate: basis set convergence check.
    energies_by_basis: dict like {'def2-svp': E1, 'def2-tzvp': E2, ...}
    Checks that the last two basis sets differ by less than threshold.
    """
    bases = list(energies_by_basis.keys())
    if len(bases) < 2:
        return True, "PASS (only one basis, no convergence check possible)"
    e_last = energies_by_basis[bases[-1]]
    e_prev = energies_by_basis[bases[-2]]
    delta_kcal = abs(e_last - e_prev) * 627.509474
    if delta_kcal > threshold_kcal:
        return False, f"BASIS_NOT_CONVERGED: Δ={delta_kcal:.2f} kcal/mol > {threshold_kcal}"
    return True, f"PASS (Δ={delta_kcal:.2f} kcal/mol)"
```

### 2.2 Gate Application Order

```
For every molecule in a pipeline:
│
├── After Step 0 (input): gate_smiles_valid
├── After Step 1 (optimization): gate_scf_converged AND gate_geometry_converged
├── After Step 2 (frequency): gate_frequency(mode="minimum")
├── After Step 3 (properties): gate_scf_converged AND gate_spin_contamination
├── After Step 4 (excited states): gate_tddft_converged (all states)
└── After Step 5 (basis convergence): gate_basis_convergence
```

## Phase 3: Cross-Agent Data Contracts

### 3.1 CE → QE: Molecule Delivery Format

When ChemicalExpert sends molecules to QuantumExpert for QC evaluation:

```json
{
  "source_agent": "ChemicalExpert",
  "task": "qc_screening",
  "project": "IPF_ALK5_cycle3",
  "molecules": [
    {
      "mol_id": "CE-IPF-001",
      "smiles": "c1ccc2c(c1)ncn2",
      "charge": 0,
      "spin": 0,
      "priority": "high",
      "ce_notes": "Top5 from docking; hinge binder confirmed"
    }
  ],
  "requested_properties": ["E_total", "HOMO", "LUMO", "gap", "dipole"],
  "method_preference": "default",
  "deadline_note": "needed for Cycle 3 decision"
}
```

**Required fields**: mol_id, smiles, charge, spin
**Optional fields**: priority, ce_notes, requested_properties, method_preference

### 3.2 QE → CE: Results Delivery Format

When QuantumExpert returns results to ChemicalExpert:

```json
{
  "source_agent": "QuantumExpert",
  "task": "qc_screening_results",
  "project": "IPF_ALK5_cycle3",
  "method_label": "PBE0-D3(BJ)/def2-TZVP, PySCF 2.12.1",
  "results": [
    {
      "mol_id": "CE-IPF-001",
      "smiles": "c1ccc2c(c1)ncn2",
      "E_total_Ha": -395.123456789,
      "HOMO_Ha": -0.2345,
      "LUMO_Ha": -0.0567,
      "gap_eV": 4.834,
      "dipole_D": 2.345,
      "ZPE_kcal": 67.89,
      "n_imag": 0,
      "qc_flag": "PASS",
      "qe_notes": "Clean convergence; no issues"
    }
  ],
  "batch_summary": {
    "n_total": 50,
    "n_passed": 47,
    "n_failed": 3,
    "failure_breakdown": {
      "SCF_FAIL": 1,
      "OPT_FAIL": 1,
      "IMAG_FREQ": 1
    }
  }
}
```

**qc_flag values**: PASS, SCF_FAIL, OPT_FAIL, IMAG_FREQ, SPIN_CONTAM, SKIPPED
**Rule**: CE decides go/no-go on drug-likeness. QE provides data + QC flags, not downstream filtering decisions.

### 3.3 QE → CE: Excited-State Results (Extended Format)

For tasks requiring electronic structure beyond ground-state properties:

```json
{
  "mol_id": "CE-OPV-042",
  "ground_state": {
    "method": "PBE0/def2-TZVP",
    "E_total_Ha": -648.229,
    "HOMO_Ha": -0.274,
    "LUMO_Ha": -0.059,
    "gap_eV": 5.85
  },
  "excited_states": {
    "method": "TD-PBE0/def2-TZVP",
    "states": [
      {
        "state": "S1",
        "E_eV": 3.42,
        "lambda_nm": 362.5,
        "f_osc": 0.8723,
        "character": "pi→pi*",
        "nto_dominant_weight": 0.94
      }
    ]
  },
  "qc_flag": "PASS",
  "qe_notes": "S1 is bright pi→pi*, good for OPV absorption"
}
```

## Phase 4: Batch Processing Strategy

### 4.1 Batch Execution Rules

```
Batch processing principles:
│
├── Independence: each molecule is processed independently
│   └── One failure does not affect others
│
├── Fail-fast per molecule: if a gate fails, skip to next molecule
│   └── Do not waste compute on downstream steps for a failed molecule
│
├── Logging: every molecule gets a result entry (PASS or specific failure code)
│   └── Never silently skip a molecule
│
├── Checkpointing: write results after EACH molecule completes
│   └── If the batch is interrupted, completed results are not lost
│   └── (Learned from CE's Cycle 1: batch docking timeout → chunk strategy)
│
├── Chunk sizing: for very large batches (>100 molecules)
│   ├── Process in chunks of 10-20
│   ├── Write intermediate batch_results_chunk_N.json after each chunk
│   └── Merge at the end
│   └── (Prevents gateway timeout; keeps memory manageable)
│
└── Progress reporting: after each chunk, print:
    └── "Chunk N/M complete: X passed, Y failed, Z remaining"
```

### 4.2 Batch Execution Script Template

```python
#!/opt/conda/envs/chem/bin/python
"""Batch QC screening pipeline template."""
import json
import time
from pathlib import Path

def process_one_molecule(mol_data, method_config, workdir):
    """Process a single molecule through the QC pipeline.

    Returns a result dict with qc_flag and all properties.
    """
    mol_id = mol_data["mol_id"]
    smiles = mol_data["smiles"]
    t0 = time.time()
    result = {"mol_id": mol_id, "smiles": smiles}

    try:
        # Step 0: Input validation
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result["qc_flag"] = "INVALID_SMILES"
            return result

        # Step 1: 3D embed
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            result["qc_flag"] = "EMBED_FAIL"
            return result
        AllChem.UFFOptimizeMolecule(mol)

        # Step 2-3: DFT optimization + properties
        # [call qchem-dft pipeline functions here]
        # ... optimization, frequency, property extraction ...

        # Gate checks
        # if not gate_scf_converged(mf): ...
        # if not gate_frequency(freqs): ...

        result["qc_flag"] = "PASS"
        result["wall_time_s"] = time.time() - t0

    except Exception as e:
        result["qc_flag"] = f"ERROR: {str(e)[:100]}"
        result["wall_time_s"] = time.time() - t0

    return result


def run_batch(input_path, output_dir, method_config, chunk_size=10):
    """Run batch QC screening with checkpointing."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(input_path) as f:
        molecules = json.load(f)["molecules"]

    all_results = []
    n_pass, n_fail = 0, 0

    for i, mol_data in enumerate(molecules):
        result = process_one_molecule(mol_data, method_config, output_dir)
        all_results.append(result)

        if result["qc_flag"] == "PASS":
            n_pass += 1
        else:
            n_fail += 1

        # Checkpoint after each molecule
        if (i + 1) % chunk_size == 0 or i == len(molecules) - 1:
            chunk_num = (i + 1) // chunk_size
            checkpoint = {
                "status": "in_progress" if i < len(molecules) - 1 else "complete",
                "processed": i + 1,
                "total": len(molecules),
                "n_pass": n_pass,
                "n_fail": n_fail,
                "results": all_results,
            }
            (output_dir / "batch_results.json").write_text(
                json.dumps(checkpoint, indent=2))
            print(f"  Checkpoint {i+1}/{len(molecules)}: "
                  f"{n_pass} passed, {n_fail} failed")

    # Final summary
    summary = {
        "n_total": len(molecules),
        "n_passed": n_pass,
        "n_failed": n_fail,
        "failure_rate": n_fail / len(molecules) if molecules else 0,
    }
    if summary["failure_rate"] > 0.2:
        print(f"  ⚠️ WARNING: failure rate {summary['failure_rate']:.0%} > 20%")
        print(f"     Something systemic may be wrong — check method/basis/input quality")

    return all_results, summary
```

## Phase 5: Audit Trail

### 5.1 Audit Log Format

Every pipeline run must produce an audit log:

```json
{
  "pipeline": "single_molecule_qc",
  "molecule": "aspirin",
  "started_at": "2026-03-03T10:00:00Z",
  "completed_at": "2026-03-03T12:30:00Z",
  "steps": [
    {
      "step": 1,
      "skill": "qchem-dft",
      "action": "geometry_optimization",
      "method": "B3LYP-D3(BJ)/def2-TZVP",
      "software": "PySCF 2.12.1",
      "gate_result": "PASS",
      "wall_time_s": 3600,
      "notes": "geomeTRIC converged in 42 steps"
    },
    {
      "step": 2,
      "skill": "qchem-dft",
      "action": "frequency_analysis",
      "method": "B3LYP-D3(BJ)/def2-TZVP",
      "software": "PySCF 2.12.1",
      "gate_result": "PASS",
      "wall_time_s": 5400,
      "notes": "n_imag=0, min_freq=16.2 cm-1"
    }
  ],
  "final_status": "COMPLETE",
  "output_files": [
    "optimized_geometry.xyz",
    "summary.json",
    "summary.md"
  ]
}
```

### 5.2 Audit Requirements

```
Every pipeline run MUST log:
├── Which skills were called, in what order
├── What method/basis/software was used at each step
├── Gate pass/fail at each step (with reason if failed)
├── Wall time per step
├── Output files produced
└── Final status (COMPLETE / PARTIAL / FAILED)

This is the lab notebook. Without it, results are not reproducible.
```

## Phase 6: Error Recovery

### 6.1 Common Interruptions and Recovery

```
Pipeline interrupted?
│
├── Codex OAuth timeout (known issue)
│   ├── Check: which step was running? Look at last checkpoint file.
│   ├── Recovery: restart from last completed step (don't redo earlier steps)
│   └── Prevention: use checkpointing; run long calculations via docker exec -d
│
├── SCF not converging for one molecule in batch
│   ├── Action: flag as SCF_FAIL, skip to next molecule
│   ├── Do NOT retry with same settings (waste of time)
│   └── If > 20% molecules fail SCF: check method/basis choice for this chemical class
│
├── Memory error during Hessian
│   ├── Action: try smaller basis for frequency (SVP instead of TZVP)
│   ├── Or: use numerical frequency with reduced displacement
│   └── Log: "freq computed at lower level due to memory constraint"
│
├── Gateway timeout on long calculation
│   ├── Prevention: run via docker exec -d (bypasses gateway)
│   ├── Recovery: check output files; restart from last checkpoint
│   └── (Same lesson as CE's docking batch: chunk to avoid timeout)
│
└── Git commit fails (no user.name/email configured)
    ├── Results are NOT lost (they're in the filesystem)
    ├── Ask user to configure git
    └── Commit after configuration is fixed
```

## Checklist Before Starting a Pipeline

- [ ] **Pipeline plan documented**: all steps, methods, gates listed before computing
- [ ] **User approved**: plan reviewed and approved (for non-trivial pipelines)
- [ ] **Input validated**: all SMILES parse, charges/spins make sense
- [ ] **Method justified**: functional/basis selection follows qchem-dft decision tree
- [ ] **Output directory created**: research/qchem/<project>/<date>/
- [ ] **Checkpointing enabled**: results saved after each step/molecule
- [ ] **Audit log initialized**: pipeline start time, parameters recorded

## Checklist After Completing a Pipeline

- [ ] **All gates checked**: every step has pass/fail logged
- [ ] **Summary written**: both JSON (machine) and MD (human)
- [ ] **Audit log complete**: all steps, timings, outcomes
- [ ] **Failures documented**: every skip/failure has a reason code
- [ ] **Results archived**: research/ directory with all artifacts
- [ ] **Git committed**: all results + scripts + logs

## Integration with Vault

- **Orchestrates**: `qchem-dft` (geometry, frequency, properties) and `qchem-excited-state` (excited states)
- **Receives from CE**: molecules via `chem-molgen` or `chem-experiment` output
- **Returns to CE**: standardized qc_results.csv/json for downstream ADMET/docking/scoring
- **Literature**: `chem-literature` for experimental reference values in benchmarking
- **Batch lessons from CE**: chunk sizing, checkpoint strategy, timeout mitigation — all inherited from CE's IPF docking experience
