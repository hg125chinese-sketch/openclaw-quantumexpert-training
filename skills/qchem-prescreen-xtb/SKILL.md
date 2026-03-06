---
name: qchem-prescreen-xtb
description: Fast pre-DFT prescreening with GFN2-xTB via the xtb binary. Use to (1) pre-optimize geometries before DFT, (2) triage VAE/generative molecules that are likely to fail SCF/optimization, (3) run optional cheap frequency checks, and (4) produce standardized prescreen QC flags and artifacts for qchem-workflow batch pipelines.
---

# xTB Prescreening (GFN2-xTB) for DFT Workflows

This skill is a **fast prescreen layer**. It does not replace DFT gates; it reduces wasted DFT time by:
- removing catastrophic geometry strain before DFT,
- providing a fail-fast heuristic for obvious problem structures,
- producing standardized artifacts + QC flags for orchestration.

## 0) Environment verification (MANDATORY gate)

**Default xTB path (this environment):** `/home/node/.local/bin/xtb`

**Gate:** xTB must be available and runnable.

Run:
```bash
/home/node/.local/bin/xtb --version || xtb --version
which xtb || true
```

If `xtb` is not found:
- **Stop** and install it (see Installation section).

## 1) Input format (data contract)

### 1.1 Single molecule JSON (recommended)

```json
{
  "mol_id": "string",
  "smiles": "optional string",
  "xyz_path": "optional path to XYZ",
  "charge": 0,
  "spin": 0,
  "settings": {
    "method": "gfn2",
    "run_freq": false,
    "max_cycles": 250,
    "etemp": 300
  }
}
```

Rules:
- Provide **either** `xyz_path` **or** `smiles`.
- `spin` uses the PySCF convention (2S). For singlet closed-shell use `spin=0`.

### 1.2 Batch JSON

```json
{
  "project": "string",
  "molecules": [ <single molecule objects> ]
}
```

## 2) Outputs (standardized)

For each molecule, write a per-molecule directory:

```
<workdir>/<mol_id>/
  input.xyz
  xtb_opt.xyz
  xtb.out
  prescreen.json
  (optional) xtb_freq.out
```

### 2.1 prescreen.json (schema)

```json
{
  "mol_id": "...",
  "method": "GFN2-xTB",
  "charge": 0,
  "spin": 0,
  "xtb_version": "...",
  "paths": {
    "input_xyz": "...",
    "xtb_out": "...",
    "opt_xyz": "..."
  },
  "timing_s": {
    "opt": 0.0,
    "freq": 0.0
  },
  "qc": {
    "embed_ok": true,
    "opt_ok": true,
    "freq_ok": false,
    "qc_flag": "PASS|XTB_EMBED_FAIL|XTB_OPT_FAIL|XTB_FREQ_FAIL",
    "notes": "free text"
  }
}
```

## 3) Procedure

### 3.1 Prepare input geometry

Preferred:
- Use an existing XYZ (from RDKit ETKDG, CE handoff, or prior workflow step).

If only SMILES is available:
- Generate a 3D geometry with RDKit ETKDG.
- Gate: embedding succeeded.

If embedding fails:
- Set `qc_flag = XTB_EMBED_FAIL` and stop processing this molecule.

### 3.2 GFN2-xTB geometry optimization (subprocess)

**Method label:** GFN2-xTB

Run (typical):
```bash
xtb input.xyz --gfn 2 --opt --chrg <charge> --uhf <uhf>
```

Notes:
- xTB uses `--uhf` as number of unpaired electrons.
  - For closed-shell singlet: `spin=0` → `uhf=0`.
  - In general: `uhf = spin` (since PySCF spin is 2S).

**Gate: OPT success**
- Process exit code == 0 (best-effort)
- AND optimized structure exists (commonly `xtbopt.xyz` in the working directory)

If optimization fails:
- Set `qc_flag = XTB_OPT_FAIL`
- Keep `xtb.out` for diagnosis
- Do not proceed to DFT.

### 3.3 Optional frequency check

If `settings.run_freq == true`, run:
```bash
xtb xtb_opt.xyz --gfn 2 --hess --chrg <charge> --uhf <uhf>
```

**Gate: FREQ success**
- Exit code == 0
- Output file exists (store as `xtb_freq.out`)

Do **not** use xTB frequencies as a replacement for the DFT `n_imag==0` gate.
Use it as a cheap diagnostic only.

If frequency fails:
- Set `qc_flag = XTB_FREQ_FAIL` (but keep `xtb_opt.xyz` as a potentially useful starting geometry for DFT if you choose to proceed).

## 4) QC flags (required)

Use exactly one of:
- `PASS` — embed + opt succeeded (freq optional)
- `XTB_EMBED_FAIL`
- `XTB_OPT_FAIL`
- `XTB_FREQ_FAIL`

## 5) Integration with qchem-workflow

Add a prescreen step **before** DFT optimization:

- Step 1: xTB prescreen
  - If `PASS`: forward `xtb_opt.xyz` as the starting geometry for `qchem-dft` / DFT geomeTRIC optimization.
  - If `XTB_OPT_FAIL`: mark as high OPT_FAIL risk; skip DFT by default.
  - If `XTB_FREQ_FAIL`: optional policy; either skip or proceed to DFT with a note.

Batch policy:
- One molecule failing does not stop the batch.
- Always checkpoint per-molecule artifacts.

## 6) Installation

Recommended (conda-forge):
```bash
conda install -c conda-forge xtb
```

Or use upstream binaries from the official releases.

## 7) Known limitations (read before using)

- **xTB PASS does not guarantee DFT PASS.** Validation on Cycle 3 (n=4) found:
  - xTB prescreen: **4/4 PASS**
  - DFT outcome: **2/4 OPT_FAIL**
  Therefore, treat xTB as a fail-fast filter and geometry preconditioner, not as a predictor of DFT convergence.

## 8) Reporting requirements

Every prescreen result must include:
- method label: **GFN2-xTB (xtb)**
- xTB version string
- qc_flag
- paths to artifacts
- wall time
