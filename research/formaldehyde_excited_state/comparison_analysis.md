# PLAYBOOK — Formaldehyde (H₂CO) Excited-State Comparison Analysis

## 0) Scope / artifacts (reproducible)

**Target:** Formaldehyde (H₂CO) — *vertical* excitations at the optimized ground-state geometry.

**Geometry / ground state method label**
- Geometry: optimized at **PBE0/def2-TZVP** (RKS), PySCF **2.12.1**
- Grid: level=4
- SCF conv_tol = 1e-10
- Geometry file: `formaldehyde_pbe0_def2-tzvp_opt.xyz`

**Excited-state methods compared (method labels)**
1) **TD-PBE0/def2-TZVP** (singlets), nstates=5, td_conv_tol=1e-6
2) **SA-CASSCF(4,3)/def2-TZVP** (6 roots = S0 + 5)
3) **NEVPT2//CASCI(4,3)/def2-TZVP @ SA-CASSCF orbitals** (per-root NEVPT2)
4) **SA-CASSCF(6,5)/def2-TZVP** (6 roots = S0 + 5; expanded active space)
5) **NEVPT2//CASCI(6,5)/def2-TZVP @ SA-CASSCF orbitals** (per-root NEVPT2)

**Primary result files**
- (4,3) + TD-DFT: `results.json`
- (6,5) + NOONs: `results_cas65.json`

**Skill evidence (active space guidance)**
From `qchem-excited-state`:
- Active space design must be validated using NOONs; for carbonyl n→π* the guide suggests **(4,3) or (6,5)** and warns it **may need σ_CO**.

```text
Source: qmd://qchem-excited-state/skill.md:246-290
268  │   ├── Natural orbital occupation numbers (NOONs):
269  │   │   ├── All occupations near 2.0 or 0.0 → orbital not needed (remove it)
270  │   │   ├── Occupations between 0.02 and 1.98 → orbital is active (keep it)
271  │   │   └── Occupations near 1.0 → strongly correlated (definitely keep)
...
287  | Carbonyl n→π* | (4,3) or (6,5) | n + π + π*; may need σ_CO |
```

---

## 1) Gates / QC status

- Ground-state SCF converged: **PASS**
- Geometry optimization converged: **PASS** (geomeTRIC “Converged!”)
- TD-DFT requested states converged: **PASS for reported S1–S5**
- SA-CASSCF(4,3) converged: **PASS**
- NEVPT2(4,3) computed for roots 0–5: **PASS**
- SA-CASSCF(6,5) converged: **PASS**
- NEVPT2(6,5) computed for roots 0–5: **PASS**

---

## 2) TD-DFT results (reference for state character)

**TD-PBE0/def2-TZVP** (PySCF 2.12.1; vertical at optimized S0 geometry):

| State | NTO type | E (eV) | λ (nm) | f |
|---:|---|---:|---:|---:|
| S1 | n→π* | 4.0429 | 306.7 | ~0.0000 |
| S2 | n→π* | 8.1786 | 151.6 | 0.1075 |
| S3 | π→π* | 9.2461 | 134.1 | 0.000171 |
| S4 | n→π* | 9.3335 | 132.8 | 0.03387 |
| S5 | n→π* (mixed; w0≈0.53) | 9.5279 | 130.1 | 0.000625 |

(NTO classification here is from transition-density SVD on the TD amplitudes.)

---

## 3) (1) TD-DFT vs NEVPT2: **state matching by excitation type** (not by index)

### 3.1 Why we match by type
The ordering of multi-reference “roots” can change when:
- the active space changes,
- PT2 corrections are state-dependent,
- some roots mix strongly.

Therefore, we match by **expected carbonyl physics**:
- The lowest singlet is typically **n→π*** around ~4 eV (experimentally).

### 3.2 Matching using the improved CAS(6,5) NEVPT2 spectrum
The **NEVPT2(6,5)** vertical excitation energies are:
- Root1–Root5: **3.6746, 4.1005, 6.3150, 10.0415, 10.2782 eV**

The first two roots (3.67–4.10 eV) now form a plausible manifold to match the experimental/TD-DFT n→π*.

**Best-effort match (by energy scale + carbonyl expectation):**
- Experimental S1 (n→π*) ≈ 4.0 eV
  - TD-DFT S1: 4.043 eV
  - NEVPT2(6,5) Root1: 3.675 eV
  - NEVPT2(6,5) Root2: 4.101 eV

At this point, the multi-reference spectrum is no longer “pathologically low” (unlike the CAS(4,3) NEVPT2 Root1 = 0.97 eV). However, *which* of Root1 vs Root2 is the n→π* vs another valence excitation still needs explicit root-character analysis (e.g., root-resolved 1-RDM differences / CI-vector diagnostics). Energetically, **Root2 (4.10 eV)** is closest to the experimental S1.

---

## 4) (2) Diagnosis of anomalous low NEVPT2 excitations: (4,3) vs (6,5)

### 4.1 PT2 correction magnitudes for CAS(4,3)
From `results.json -> multiref.nevpt2.corr_ha`:
- corr(root0..5) = [-0.4249, -0.4602, -0.3993, -0.4779, -0.4133, -0.4344] Ha

Key diagnostic (state-dependent PT2):
- Δcorr(root1 − root0) ≈ **-0.0353 Ha ≈ -0.96 eV**
- Δcorr(root3 − root0) ≈ **-0.0530 Ha ≈ -1.44 eV**

This large *differential* PT2 stabilization reorders roots and produced an unphysical-looking very low excitation (Root1 ≈ 0.97 eV).

### 4.2 PT2 correction magnitudes for CAS(6,5)
From `results_cas65.json -> nevpt2.corr_ha`:
- corr(root0..5) = [-0.36867, -0.37449, -0.37418, -0.35510, -0.38143, -0.38188] Ha

Convert spread:
- max(corr) − min(corr) ≈ (-0.35510) − (-0.38188) = **0.02678 Ha ≈ 0.73 eV**

This is **much more uniform** than the CAS(4,3) case (spread ≈ 0.07865 Ha ≈ 2.14 eV). The expanded active space materially reduces the root-to-root variability in PT2 stabilization.

### 4.3 NOON validation (CAS(6,5))
From `results_cas65.json -> casci.noons` (per root, sorted descending):
- Root0 NOONs: [1.999, 1.981, 1.925, 0.076, 0.019]
  - Interpretable as: 3 mostly doubly-occupied active orbitals and 2 mostly virtual → consistent with a closed-shell reference within the chosen valence space.
- Root1 NOONs: [1.988, 1.982, 1.011, 1.008, 0.010]
- Root2 NOONs: [1.988, 1.981, 1.011, 1.010, 0.010]
- Root3 NOONs: [1.999, 1.991, 1.000, 1.000, 0.010]
  - Roots 1–3 show two orbitals near ~2, two near ~1, one near ~0 → consistent with “one-electron promotion” physics in the active space.
- Roots 4–5 NOONs: ~[1.85, 1.83, 1.15, 1.15, 0.02]
  - More strongly mixed/correlated within the active space (occupations closer to 1), consistent with higher excited states and/or stronger configuration mixing.

**Conclusion of NOON gate:** the (6,5) active space is *active* (has nontrivial occupations away from 0/2, especially for excited roots), and therefore is more likely to provide a stable NEVPT2 reference than (4,3).

---

## 5) (3) Experimental check for S1 (n→π*) ≈ 4.0 eV

User reference: **S1 ≈ 4.0 eV (n→π*)**.

- TD-PBE0/def2-TZVP S1 = **4.0429 eV**
  - Deviation ≈ **+0.043 eV** (excellent vs the provided reference)

With the improved CAS(6,5):
- NEVPT2 Root2 = **4.1005 eV** (Δ ≈ +0.10 eV)
- NEVPT2 Root1 = **3.6746 eV** (Δ ≈ -0.33 eV)

Given only energy information, **Root2** is the most plausible match to the experimental n→π* S1. To make this definitive, we should compute a root-character diagnostic (e.g., difference density / orbital transition analysis) for the CASCI roots.

---

## 6) Bottom line / what changed after expanding to (6,5)

1) The suspiciously low **NEVPT2(4,3) Root1 = 0.97 eV** disappears. The lowest NEVPT2(6,5) excitations are now **3.67–4.10 eV**, consistent with carbonyl valence physics.
2) PT2 corrections become **more uniform** across roots (spread drops from ~2.14 eV equivalent to ~0.73 eV equivalent), indicating the reference space is less “imbalanced” between states.
3) Root-resolved NOONs provide internal validation that the (6,5) space is actually capturing active correlation and excitation character.

---

## 7) Next optional diagnostics (if you want a stricter state-by-type mapping)

- Compute a per-root difference density (CASCI root − ground) in AO basis and quantify charge transfer / lone pair depletion on O.
- Alternatively, compute natural orbitals per root and identify which orbital corresponds to n(O) vs π/σ.

That would let us label NEVPT2 roots explicitly as n→π* vs π→π* and match to the TD-DFT states without ambiguity.
