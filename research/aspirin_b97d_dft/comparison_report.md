# PLAYBOOK — Aspirin DFT Comparison Report (PBE0/no-disp vs B97-D/disp-inclusive)

## 0) Scope / artifacts (reproducible)

**Molecule**: aspirin (acetylsalicylic acid)  
SMILES: `CC(=O)Oc1ccccc1C(=O)O`  
Charge = 0; Spin = 0 (RKS)

**Software / numerical settings**
- PySCF: 2.12.1
- libxc: 7.0.0 (via PySCF)
- Grid: level = 4 (Treutler–Ahlrichs radial + Becke partition; as printed by PySCF)
- SCF conv_tol = 1e-9

**Compared runs (method labels)**
1) **PBE0/def2-SVP → def2-TZVP**, *no dispersion* (tooling-limited; D3 unavailable)  
   Source: `../aspirin_dft/summary.json`
2) **B97-D/def2-SVP → def2-TZVP**, *dispersion-inclusive XC* (built-in dispersion in the XC definition)  
   Source: `summary.json`

**Important comparability note (PLAYBOOK)**
- Absolute energies across different basis sets and/or separately optimized geometries are **not** used as a basis-set convergence conclusion.
- For basis convergence, we emphasize **property stability** (gap, dipole, ZPE, and frequency gate outcomes).

---

## 1) Gates / QC status

All cases below passed the workflow gates:
- SCF converged: **true**
- Geometry optimization: **Converged** (optimized XYZ written)
- Frequency gate: **n_imag = 0**

---

## 2) Key results (method-labeled)

### 2.1 HOMO–LUMO gap (KS gap)

| Method label | HOMO (Ha) | LUMO (Ha) | Gap (eV) |
|---|---:|---:|---:|
| **PBE0/def2-SVP** (PySCF 2.12.1; no-disp) | -0.27040180 | -0.05557668 | **5.8457** |
| **PBE0/def2-TZVP** (PySCF 2.12.1; no-disp) | -0.27428441 | -0.05923097 | **5.8519** |
| **B97-D/def2-SVP** (PySCF 2.12.1; built-in disp) | -0.22440909 | -0.08501353 | **3.7931** |
| **B97-D/def2-TZVP** (PySCF 2.12.1; built-in disp) | -0.23080538 | -0.08997935 | **3.8321** |

### 2.2 Dipole moment (norm)

| Method label | Dipole (Debye) |
|---|---:|
| **PBE0/def2-SVP** (no-disp) | 1.8885 |
| **PBE0/def2-TZVP** (no-disp) | 1.9942 |
| **B97-D/def2-SVP** (built-in disp) | 1.8520 |
| **B97-D/def2-TZVP** (built-in disp) | 1.9912 |

### 2.3 ZPE and minimum frequency

| Method label | n_imag | min freq (cm⁻¹) | ZPE (kcal/mol) |
|---|---:|---:|---:|
| **PBE0/def2-SVP** (no-disp) | 0 | 35.83 | 99.72 |
| **PBE0/def2-TZVP** (no-disp) | 0 | 24.65 | 99.02 |
| **B97-D/def2-SVP** (built-in disp) | 0 | 26.84 | 96.36 |
| **B97-D/def2-TZVP** (built-in disp) | 0 | 16.17 | 95.74 |

---

## 3) Analysis focus

### (1) Source of the gap differences (functional vs dispersion vs basis)

**Observed**: PBE0 gaps (~5.85 eV) are ~**2.0 eV larger** than B97-D gaps (~3.8 eV).

**Primary driver: the XC functional (hybrid vs GGA)**
- **PBE0** is a **hybrid** functional (includes a fixed fraction of exact exchange). This typically **increases** the KS HOMO–LUMO gap relative to GGAs.
- **B97-D** (as used here via libxc in PySCF) is **GGA-class** with an empirical dispersion correction embedded in the functional definition. GGAs commonly **underestimate** KS gaps due to self-interaction/delocalization error and the lack of nonlocal exchange.
- The magnitude here (~2 eV) is consistent with a *functional-class* effect and is therefore attributed mainly to **functional choice**.

**Dispersion contribution: usually second-order for KS gaps**
- Dispersion corrections primarily affect **intermolecular / intramolecular weak interactions** and thus can shift **optimized geometries**.
- For a closed-shell organic molecule, the direct impact of dispersion on **frontier orbital energies** (and thus KS gap) is generally modest compared to switching between hybrid and GGA.
- In this comparison, dispersion strategy is *confounded* with functional choice (PBE0/no-disp vs B97-D/built-in-disp), but the scale of the gap change strongly points to the **functional**, not dispersion, as the dominant source.

**Basis contribution: small within each functional**
- PBE0: gap changes by **+0.0062 eV** (SVP → TZVP)
- B97-D: gap changes by **+0.0389 eV** (SVP → TZVP)
- Conclusion: for this system, **basis effects on the KS gap are minor** relative to functional effects.

### (2) Basis-set convergence assessment (def2-SVP → def2-TZVP)

**Gap stability**
- Both functionals show **very small** SVP→TZVP gap changes (≤0.04 eV). For reporting a KS gap at this level, **def2-TZVP** appears effectively converged with respect to this property.

**Dipole stability**
- PBE0: Δμ ≈ **+0.106 D**
- B97-D: Δμ ≈ **+0.139 D**
- This is a noticeable but not extreme basis effect; **TZVP** is recommended when the dipole is a key deliverable.

**Vibrational / ZPE stability**
- PBE0: ZPE decreases by **~0.71 kcal/mol**
- B97-D: ZPE decreases by **~0.62 kcal/mol**
- This indicates mild basis sensitivity in harmonic frequencies/ZPE (often influenced by low-frequency modes). TZVP is more robust.

**Total energies**
- Large decreases in total energy from SVP→TZVP are expected, but per PLAYBOOK we **do not** use absolute energy differences between separately optimized geometries as a convergence claim.

### (3) Method-selection recommendations for future similar tasks

**If the deliverable emphasizes an interpretable/reproducible gap**
- Avoid using a **GGA** (including B97-D) as the final “gap answer”. GGAs typically underestimate KS gaps.
- Prefer:
  - **Hybrid** (e.g., PBE0) for a stable, commonly used KS-gap baseline
  - Or **range-separated hybrids** when available/desired for improved frontier orbital energetics (requires checking PySCF/libxc support in the environment)

**Dispersion compliance under tooling limits (D3/D4 unavailable)**
- Continue following the rule: *do not silently drop dispersion*.
- If external D3 cannot be used, a compliant fallback is to switch to an **XC that includes dispersion/long-range correlation** in its definition (e.g., B97-D, or other supported options).

**Practical split-workflow recommendation (reduces confounding)**
- Use a **dispersion-inclusive method** for geometry + frequencies (structure-sensitive steps).
- Then compute **gap-relevant single points** on the *same optimized geometry* with a **hybrid / RS-hybrid** (minimizes geometry and basis confounding when comparing gaps).

**Basis recommendation**
- **def2-SVP**: good for fast screening / initial optimization.
- **def2-TZVP**: recommended default for final reporting of dipoles, frequencies/ZPE, and for more stable comparisons.

---

## 4) Bottom line (executive summary)

1) The ~2 eV gap difference between PBE0 and B97-D is dominated by **functional choice (hybrid vs GGA)**; basis effects are small and dispersion is expected to be secondary for KS gap.  
2) def2-SVP → def2-TZVP shows **good practical convergence** for the KS gap (≤0.04 eV shift) and modest improvements for dipole/ZPE.  
3) For future work: keep **dispersion on** for geometries/frequencies (via available tooling or dispersion-inclusive XC), and use **hybrid/RS-hybrid single points** for gap-focused deliverables.
