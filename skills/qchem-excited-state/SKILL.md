---
name: qchem-excited-state
description: Excited-state calculations including TD-DFT, CASSCF, NEVPT2, and CASPT2. Covers method selection, active space design, state-averaging, solvatochromism, and the critical decision of when single-reference methods fail and multi-reference is required.
homepage: https://pyscf.org
metadata: { "openclaw": { "emoji": "✨", "requires": { "bins": ["python3"], "python": ["pyscf", "numpy", "scipy", "matplotlib"] } } }
---

# Excited-State Calculations: TD-DFT & Multi-Reference Methods

Covers the full spectrum of excited-state quantum chemistry: from fast TD-DFT screening to rigorous multi-reference (CASSCF/NEVPT2/CASPT2) treatments. The central challenge is knowing **when single-reference fails** and **how to build the right active space** — this skill provides decision frameworks for both.

## When to Use

- User asks for UV-Vis absorption / emission wavelengths or oscillator strengths
- User needs excited-state geometry optimization or emission energy
- User asks about charge-transfer (CT) states, π→π*, n→π* transitions
- User wants to compare TD-DFT vs multi-reference results
- User needs to design an active space for CASSCF/NEVPT2
- User asks about conical intersections or photochemistry
- User has a system where TD-DFT is expected to fail (CT, double excitation, diradical)
- Upstream from qchem-dft: excited-state calculations require a converged ground-state first

## Core Philosophy

1. **TD-DFT first, multi-reference when needed.** TD-DFT is cheap and often good enough. But know its failure modes — CT states, double excitations, and near-degeneracies will silently give wrong answers.
2. **Active space design is the bottleneck.** CASSCF is only as good as your active space. Too small → missing physics. Too large → intractable. This is where chemical intuition meets computational reality.
3. **State-averaging is not optional for CASSCF.** If you want excited states, you must state-average. Single-root CASSCF for excited states is wrong — it biases toward the ground state.
4. **Vertical ≠ adiabatic ≠ 0-0.** Always specify which energy you're computing and which one the experiment measures.
5. **Benchmark against experiment with the right quantity.** λ_max from UV-Vis corresponds to vertical excitation. Fluorescence λ_em corresponds to vertical emission from excited-state geometry. 0-0 transitions need both geometries + ZPE correction.

## Phase 1: Method Selection Decision Tree

```
What excited-state calculation do you need?
│
├── Quick screening / large system (> 50 atoms)
│   └── TD-DFT (Phase 2)
│       ├── Local excitation (π→π*, n→π*) → works well
│       ├── Charge-transfer state → use range-separated functional (ωB97X-D, CAM-B3LYP)
│       └── Rydberg state → add diffuse functions (aug-cc-pVTZ or ma-def2-TZVP)
│
├── Medium system (15–50 atoms), need reliable energetics
│   ├── Is it single-reference? (check Λ diagnostic, see Phase 2.4)
│   │   ├── YES (Λ > 0.3) → TD-DFT is likely fine
│   │   └── NO (Λ < 0.3) or unsure → consider multi-reference
│   └── Does it involve double excitation character?
│       ├── NO → TD-DFT may still work
│       └── YES → TD-DFT will fail → multi-reference required (Phase 3)
│
├── Small system (< 15 atoms), need high accuracy
│   ├── Single-reference character → EOM-CCSD (if available) or ADC(2)
│   └── Multi-reference character → CASSCF + NEVPT2 (Phase 3)
│
├── Photochemistry / conical intersections
│   └── CASSCF mandatory (TD-DFT cannot describe CI topology)
│
└── Near-degenerate ground state / diradical
    └── CASSCF + NEVPT2 mandatory (Phase 3)
```

## Phase 2: TD-DFT

### 2.1 Standard TD-DFT Calculation

```python
#!/usr/bin/env python
"""TD-DFT vertical excitation energies with PySCF."""
from pyscf import gto, dft, tddft
import numpy as np

mol = gto.M(
    atom='''...''',  # optimized ground-state geometry (from qchem-dft)
    basis='def2-tzvp',
    charge=0,
    spin=0,
    verbose=4,
)

# Ground-state SCF (must converge first — gate from qchem-dft)
mf = dft.RKS(mol)
mf.xc = 'pbe0'  # or cam-b3lyp for CT states, wb97x-d for general use
mf.grids.level = 4
mf.conv_tol = 1e-9
mf.kernel()
assert mf.converged, "Ground-state SCF not converged — cannot proceed to TD-DFT"

# TD-DFT
td = tddft.TDDFT(mf)
td.nstates = 10        # number of excited states
td.conv_tol = 1e-6     # convergence for Davidson solver
td.kernel()

# Extract results
print("\n=== TD-DFT Excitation Energies ===")
print(f"Method: TD-{mf.xc}/def2-tzvp | Software: PySCF {gto.__version__}")
print(f"{'State':>6} {'E (eV)':>10} {'λ (nm)':>10} {'f':>10} {'Character':>15}")
print("-" * 55)

for i, (e_ev, f_osc) in enumerate(zip(td.e * 27.211386, td.oscillator_strength())):
    lam_nm = 1239.84 / e_ev if e_ev > 0 else float('inf')
    print(f"  S{i+1:>3d}  {e_ev:10.4f}  {lam_nm:10.1f}  {f_osc:10.4f}")
```

### 2.2 Functional Selection for Excited States

```
Excited-state functional selection:
│
├── Local valence excitation (π→π*, n→π*)
│   ├── B3LYP — decent, systematic errors ~0.2–0.3 eV
│   ├── PBE0 — similar quality, slightly better for some systems
│   └── Both are acceptable starting points
│
├── Charge-transfer excitation
│   ├── B3LYP / PBE0 → WILL FAIL (underestimate CT energy by 1–2 eV)
│   ├── CAM-B3LYP — range-separated, much better for CT
│   ├── ωB97X-D — good all-rounder with RS correction
│   └── Rule: if donor-acceptor distance > 3 Å, use range-separated
│
├── Rydberg states
│   └── Any functional + diffuse basis (aug-cc-pVTZ or ma-def2-TZVP)
│       Standard def2-TZVP lacks diffuse functions → Rydberg states artifactually high
│
└── If unsure
    └── ωB97X-D or CAM-B3LYP (range-separated handles most cases reasonably)
```

### 2.3 Basis Set Considerations for Excited States

```
Basis set for excited states:
│
├── Valence excitations → def2-TZVP (usually sufficient)
│
├── Rydberg / diffuse states → MUST add diffuse functions
│   ├── aug-cc-pVTZ (Dunning augmented)
│   ├── ma-def2-TZVP (minimally augmented — cheaper, often enough)
│   └── Without diffuse: Rydberg states appear at wrong energy or missing entirely
│
├── Charge-transfer → def2-TZVP usually OK (CT is valence-like)
│
└── Basis set convergence for excited states:
    - Excitation energies converge faster than total energies
    - TZ is usually sufficient; QZ rarely needed for TD-DFT
    - Check: does adding diffuse functions change your state ordering?
```

### 2.4 TD-DFT Diagnostic: When to Worry

```python
def tddft_diagnostics(td, threshold_lambda=0.3):
    """Check TD-DFT reliability diagnostics.

    Lambda diagnostic (Peach et al. 2008):
    - Λ > 0.3: likely reliable (local excitation)
    - Λ < 0.3: possible CT problem with non-RS functional
    - Λ < 0.15: almost certainly unreliable without RS functional

    Also check for:
    - Very low oscillator strengths (dark states — may need more states)
    - States with large doubles character (TD-DFT cannot capture)
    """
    print("\n=== TD-DFT Diagnostics ===")

    for i, e_ev in enumerate(td.e * 27.211386):
        f_osc = td.oscillator_strength()[i]

        # Approximate Λ from transition density
        # (Full Λ requires NTO analysis; this is a simplified check)
        # For now, flag based on excitation character
        print(f"  S{i+1}: E={e_ev:.3f} eV, f={f_osc:.4f}")

        if f_osc < 1e-4:
            print(f"    ⚠️ Dark state (f ≈ 0) — check symmetry or state character")

    print(f"\n  Recommendation:")
    print(f"  - If using B3LYP/PBE0 and CT states suspected: switch to CAM-B3LYP/ωB97X-D")
    print(f"  - If states show strong doubles character: TD-DFT unreliable → use CASSCF/NEVPT2")
    print(f"  - Compare nstates=10 vs nstates=20 to check if state ordering is stable")
```

### 2.5 Natural Transition Orbitals (NTOs)

```python
def compute_ntos(td, state_idx=0):
    """Compute Natural Transition Orbitals for excited state analysis.

    NTOs give the most compact description of the excitation:
    hole → particle orbital pairs, ranked by singular value.
    """
    from pyscf.tdscf import rhf as td_rhf

    # Get transition density matrix
    # PySCF provides NTO analysis
    weights, nto_coeff = td.get_nto(state=state_idx + 1)

    print(f"\n=== NTO Analysis for S{state_idx+1} ===")
    print(f"{'Pair':>6} {'Weight':>10} {'Character':>20}")
    print("-" * 40)
    for i, w in enumerate(weights):
        if w > 0.05:  # only show significant pairs
            print(f"  {i+1:>4d}  {w:10.4f}")

    dominant_weight = weights[0] if len(weights) > 0 else 0
    if dominant_weight > 0.9:
        print(f"  → Single-excitation dominated (weight {dominant_weight:.3f})")
        print(f"  → TD-DFT likely reliable for this state")
    elif dominant_weight > 0.7:
        print(f"  → Mostly single-excitation (weight {dominant_weight:.3f})")
        print(f"  → TD-DFT probably OK but verify with multi-ref if critical")
    else:
        print(f"  → Significant multi-excitation character (weight {dominant_weight:.3f})")
        print(f"  → TD-DFT UNRELIABLE — use CASSCF/NEVPT2")

    return weights, nto_coeff
```

## Phase 3: Multi-Reference Methods (CASSCF / NEVPT2)

### 3.1 When Multi-Reference is Required

```
Is multi-reference needed?
│
├── Definitely YES:
│   ├── Bond breaking / dissociation
│   ├── Conical intersections / photochemistry
│   ├── Diradical / biradical character
│   ├── Near-degenerate states (e.g., transition metal open-shell)
│   ├── Double excitations (TD-DFT cannot describe)
│   └── Λ diagnostic < 0.15 with non-RS functional
│
├── Probably YES:
│   ├── Excited states of extended π-systems with low-lying dark states
│   ├── Metal-to-ligand charge transfer (MLCT) in TM complexes
│   └── Systems where TD-DFT gives qualitatively different state ordering
│       with different functionals
│
└── Probably NO (TD-DFT likely sufficient):
    ├── Simple π→π* in small aromatics
    ├── n→π* in carbonyls / azines
    ├── Single-excitation dominated (NTO weight > 0.9)
    └── Λ > 0.3 with hybrid functional
```

### 3.2 Active Space Design

This is the **most critical step** in any multi-reference calculation. A bad active space gives meaningless results.

```
Active space design strategy:
│
├── Step 1: Identify the physics
│   ├── What bonds break? → include bonding + antibonding orbitals
│   ├── What excitations matter? → include occupied + virtual orbitals involved
│   ├── Lone pairs participating? → include them
│   └── π-system? → include all π + π* orbitals
│
├── Step 2: Start small, grow carefully
│   ├── Minimum viable: (2,2) — 2 electrons in 2 orbitals
│   ├── Typical organic π→π*: (4,4) to (10,10)
│   ├── Transition metal d-shell: (n_d_electrons, 5) minimum
│   │   └── Better: include double-d shell → (n_d, 10)
│   ├── DO NOT exceed (16,16) without DMRG
│   └── Rule of thumb: N_det ∝ C(2*n_orb, n_elec) — check this is tractable
│
├── Step 3: Validate the active space
│   ├── Natural orbital occupation numbers (NOONs):
│   │   ├── All occupations near 2.0 or 0.0 → orbital not needed (remove it)
│   │   ├── Occupations between 0.02 and 1.98 → orbital is active (keep it)
│   │   └── Occupations near 1.0 → strongly correlated (definitely keep)
│   ├── Run CASSCF → check NOONs → remove dead orbitals → re-run
│   └── Energy should be stable when adding/removing borderline orbitals
│
└── Step 4: Document your choice
    └── Always report: (n_elec, n_orb), which orbitals, why, NOONs
```

### 3.3 Active Space Size Guide

| System type | Recommended CAS | Notes |
|------------|----------------|-------|
| Ethylene π→π* | (2,2) | Minimal; textbook example |
| Butadiene | (4,4) | Full π-space |
| Benzene | (6,6) | Full π-space |
| Naphthalene | (10,10) | Full π-space; getting expensive |
| Carbonyl n→π* | (4,3) or (6,5) | n + π + π*; may need σ_CO |
| Fe(II) d⁶ complex | (6,5) min; (6,10) better | d-shell; double-d for quantitative |
| Cr₂ (bond breaking) | (12,12) | Classic multi-ref benchmark |
| > (14,14) | Consider DMRG | Conventional CASSCF intractable |

### 3.4 CASSCF in PySCF

```python
#!/usr/bin/env python
"""State-averaged CASSCF for excited states with PySCF."""
from pyscf import gto, scf, mcscf
import numpy as np

mol = gto.M(
    atom='''...''',  # optimized geometry
    basis='def2-tzvp',
    charge=0,
    spin=0,
    verbose=4,
)

# Step 1: Ground-state HF (starting orbitals for CASSCF)
mf = scf.RHF(mol)
mf.kernel()
assert mf.converged, "HF not converged — cannot start CASSCF"

# Step 2: Define active space
n_elec = 6    # electrons in active space
n_orb = 6     # orbitals in active space
# Example: benzene (6,6) — all π and π* orbitals

mc = mcscf.CASSCF(mf, n_orb, n_elec)

# Step 3: State-averaging (MANDATORY for excited states)
# Average over ground state + N excited states with equal weights
n_states = 4  # S0 + 3 excited states
weights = [1.0 / n_states] * n_states
mc = mc.state_average_(weights)

# Step 4: Run CASSCF
mc.kernel()

if not mc.converged:
    print("⚠️ CASSCF not converged — try:")
    print("  1) Different initial orbitals (DFT natural orbitals)")
    print("  2) Larger max_cycle")
    print("  3) Different active space")
else:
    print("✅ SA-CASSCF converged")

# Step 5: Analyze results
print(f"\n=== SA-CASSCF({n_elec},{n_orb}) / {mol.basis} ===")
print(f"States averaged: {n_states}, weights: {weights}")
print(f"Software: PySCF {gto.__version__}\n")

for i, e in enumerate(mc.e_states):
    de = (e - mc.e_states[0]) * 27.211386  # eV relative to S0
    print(f"  State {i}: E = {e:.10f} Ha  ΔE = {de:.4f} eV")

# Step 6: Check natural orbital occupation numbers (NOONs)
# Critical for validating active space
from pyscf.tools import molden
natocc = mc.mo_occ  # only meaningful for single-state; for SA, use:

for i in range(n_states):
    dm = mc.fcisolver.make_rdm1(mc.ci[i], n_orb, n_elec)
    noons = np.sort(np.linalg.eigvalsh(dm))[::-1]
    print(f"\n  State {i} NOONs: {np.array2string(noons, precision=4)}")
    dead = sum(1 for n in noons if n > 1.98 or n < 0.02)
    if dead > 0:
        print(f"    ⚠️ {dead} orbital(s) with occupation near 0 or 2 — consider removing")
```

### 3.5 Improving CASSCF: Initial Orbital Selection

```
CASSCF not converging or giving wrong states?
│
├── Problem: HF orbitals are a bad starting point
│   └── Fix: Use DFT natural orbitals instead
│
│   # Generate DFT natural orbitals as CASSCF starting point:
│   mf_dft = dft.RKS(mol)
│   mf_dft.xc = 'pbe0'
│   mf_dft.kernel()
│   # Use DFT orbitals to initialize CASSCF:
│   mo_init = mcscf.project_init_guess(mc, mf_dft.mo_coeff)
│   mc.kernel(mo_init)
│
├── Problem: Active space orbitals rotate to wrong space
│   └── Fix: Use AVAS (Atomic Valence Active Space) for automated selection
│
│   from pyscf import mcscf
│   # AVAS: project onto specific AO character (e.g., p-orbitals for π-system)
│   avas_obj = mcscf.avas.AVAS(mf, ['C 2p', 'N 2p', 'O 2p'])
│   n_orb, n_elec, mo_coeff = avas_obj.kernel()
│   mc = mcscf.CASSCF(mf, n_orb, n_elec)
│   mc.kernel(mo_coeff)
│
└── Problem: States swap ordering during optimization
    └── Fix: Increase n_states in state-averaging; check with different weights
```

### 3.6 NEVPT2: Dynamic Correlation on Top of CASSCF

```python
#!/usr/bin/env python
"""NEVPT2 on top of SA-CASSCF for quantitative excited-state energies."""
from pyscf import gto, scf, mcscf, mrpt

mol = gto.M(atom='''...''', basis='def2-tzvp', charge=0, spin=0, verbose=4)
mf = scf.RHF(mol)
mf.kernel()

# SA-CASSCF first
n_elec, n_orb, n_states = 6, 6, 4
mc = mcscf.CASSCF(mf, n_orb, n_elec)
mc = mc.state_average_([1.0/n_states] * n_states)
mc.kernel()
assert mc.converged, "SA-CASSCF must converge before NEVPT2"

# NEVPT2 for each state
print(f"\n=== NEVPT2({n_elec},{n_orb}) / {mol.basis} ===")
print(f"Method: SC-NEVPT2 on SA-CASSCF({n_elec},{n_orb})")
print(f"Software: PySCF {gto.__version__}\n")

e_nevpt2 = []
for i in range(n_states):
    e_corr = mrpt.NEVPT(mc, root=i).kernel()
    e_total = mc.e_states[i] + e_corr
    e_nevpt2.append(e_total)
    de_cas = (mc.e_states[i] - mc.e_states[0]) * 27.211386
    de_nev = (e_total - e_nevpt2[0]) * 27.211386 if i > 0 else 0.0
    print(f"  State {i}: CASSCF ΔE = {de_cas:.4f} eV | NEVPT2 ΔE = {de_nev:.4f} eV | PT2 corr = {e_corr:.6f} Ha")

print(f"\nNEVPT2 shifts excitation energies by including dynamic correlation.")
print(f"Typical effect: 0.1–0.5 eV shift relative to CASSCF.")
print(f"If shift > 1.0 eV → active space may be too small (missing important correlations).")
```

### 3.7 NEVPT2 vs CASPT2: Which to Use?

```
NEVPT2 vs CASPT2:
│
├── NEVPT2 (recommended default in PySCF)
│   ├── No intruder state problem (numerically stable)
│   ├── Available in PySCF natively
│   ├── Slightly less accurate than CASPT2 for some benchmarks
│   └── Use this unless you have a specific reason for CASPT2
│
├── CASPT2
│   ├── More established in literature (older method)
│   ├── Can suffer from intruder states → needs IPEA shift (0.25 default)
│   ├── NOT natively available in PySCF (need OpenMolcas or BAGEL)
│   └── Use if: reproducing literature CASPT2 results, or NEVPT2 gives suspicious energies
│
└── Decision:
    ├── PySCF environment → NEVPT2 (no choice needed)
    ├── ORCA available → both available; start with NEVPT2
    └── OpenMolcas available → CASPT2 with IPEA=0.25
```

## Phase 4: Solvatochromism & Environmental Effects

### 4.1 Solvent Effects on Excited States

```python
def tddft_with_solvent(mol, xc='cam-b3lyp', basis='def2-tzvp',
                        nstates=10, eps=78.39):
    """TD-DFT with PCM solvation for solvatochromic shifts.

    Important: Linear-response PCM (LR-PCM) is the standard for TD-DFT.
    State-specific PCM (SS-PCM) is more accurate but expensive.
    """
    from pyscf import dft, tddft, solvent

    mf = dft.RKS(mol)
    mf.xc = xc
    mf = solvent.ddCOSMO(mf)
    mf.with_solvent.eps = eps
    mf.kernel()
    assert mf.converged

    td = tddft.TDDFT(mf)
    td.nstates = nstates
    td.kernel()

    print(f"\n=== TD-{xc}/{basis} + ddCOSMO(ε={eps}) ===")
    for i, (e_ev, f) in enumerate(zip(td.e * 27.211386, td.oscillator_strength())):
        lam = 1239.84 / e_ev if e_ev > 0 else float('inf')
        print(f"  S{i+1}: {e_ev:.4f} eV ({lam:.1f} nm), f={f:.4f}")

    return td
```

### 4.2 Comparing Gas Phase vs Solvent

```
Solvatochromism analysis protocol:
1) Run TD-DFT in gas phase → record excitation energies
2) Run TD-DFT with PCM (target solvent) → record excitation energies
3) Δ = E_solvent - E_gas
   - Δ < 0 (red shift / bathochromic): state stabilized by solvent
   - Δ > 0 (blue shift / hypsochromic): state destabilized by solvent
4) Compare to experimental solvatochromic shift
5) If |Δ_calc - Δ_exp| > 0.3 eV → check:
   - Is explicit solvent needed? (H-bonding solvents)
   - Is state-specific PCM needed? (large density change on excitation)
```

## Phase 5: Excited-State Geometry Optimization

### 5.1 TD-DFT Excited-State Optimization

```python
#!/usr/bin/env python
"""Optimize geometry on excited-state surface (TD-DFT)."""
from pyscf import gto, dft, tddft
from pyscf.geomopt.geometric_solver import optimize

mol = gto.M(atom='''...''', basis='def2-tzvp', charge=0, spin=0)

mf = dft.RKS(mol)
mf.xc = 'pbe0'
mf.kernel()

td = tddft.TDDFT(mf)
td.nstates = 5

# Optimize S1 geometry
td_scanner = td.as_scanner()
td_scanner.state = 1  # S1 (0-indexed in some interfaces; check PySCF docs)

mol_s1 = optimize(td_scanner, maxsteps=100)

# Vertical emission: S1 energy at S1 geometry
td_s1 = tddft.TDDFT(dft.RKS(mol_s1).run())
td_s1.nstates = 5
td_s1.kernel()

e_emission = td_s1.e[0] * 27.211386  # S1→S0 at S1 geometry
lam_em = 1239.84 / e_emission
print(f"Vertical emission: {e_emission:.4f} eV ({lam_em:.1f} nm)")
```

### 5.2 Stokes Shift

```
Stokes shift = E_absorption - E_emission

Protocol:
1) Optimize S0 geometry → compute vertical absorption (S0 geom, S1 energy)
2) Optimize S1 geometry → compute vertical emission (S1 geom, S1→S0 energy)
3) Stokes shift = E_abs - E_em

Typical Stokes shifts:
- Rigid aromatics: 0.1–0.3 eV (small geometry change)
- Flexible molecules: 0.3–1.0 eV (large relaxation)
- ESIPT (excited-state intramolecular proton transfer): > 1.0 eV
```

## Phase 6: Resource Assessment & Practical Limits

### 6.1 Feasibility Check Before Computing

```
Before starting an excited-state calculation, estimate:
│
├── TD-DFT:
│   ├── Cost ≈ 3–5× ground-state DFT per state
│   ├── Memory: similar to ground state (iterative Davidson)
│   ├── 50 atoms / def2-TZVP / 10 states → minutes to hours on laptop
│   └── 200+ atoms → still feasible with RI approximation
│
├── CASSCF:
│   ├── Cost scales as O(n_orb^6 × n_det) where n_det = C(2*n_orb, n_elec)
│   ├── (6,6): n_det = 400 → seconds
│   ├── (10,10): n_det ~63,000 → minutes
│   ├── (14,14): n_det ~40 million → hours, needs good machine
│   ├── (16,16): n_det ~600 million → may need days or DMRG
│   └── > (16,16): conventional CASSCF intractable → DMRG-CASSCF required
│
├── NEVPT2 on top of CASSCF:
│   ├── Additional cost: O(n_basis^4 × n_det) per state
│   ├── Memory: can be significant for large basis + large CAS
│   └── (10,10)/def2-TZVP on 20-atom molecule: ~1–4 hours, ~16–32 GB RAM
│
└── Memory rules of thumb:
    ├── 8 GB RAM: TD-DFT up to ~100 atoms; CASSCF up to (10,10)
    ├── 16 GB RAM: TD-DFT up to ~200 atoms; CASSCF up to (12,12)
    ├── 64 GB RAM: CASSCF (14,14); NEVPT2 on (12,12) with TZ basis
    └── > 128 GB: needed for (16,16)+ or large basis NEVPT2
```

### 6.2 When Resources Are Insufficient

```
Can't afford the "right" method?
│
├── Step 1: Can you shrink the active space?
│   └── Remove orbitals with NOONs > 1.98 or < 0.02
│
├── Step 2: Can you use a smaller basis for CASSCF?
│   └── Optimize active space at SVP, then NEVPT2 at TZ on same orbitals
│
├── Step 3: Can you truncate the molecule?
│   └── Replace distant substituents with H; keep chromophore intact
│
├── Step 4: Fallback to TD-DFT with diagnostics
│   └── Use range-separated functional + NTO analysis
│   └── Report as "TD-DFT (single-reference diagnostic: Λ = X.XX)"
│
└── Step 5: Honestly report limitations
    └── "CASPT2/NEVPT2 recommended but infeasible with current resources.
         TD-DFT results reported with caveats."
    DO NOT silently use an inadequate method.
    (Same principle as dispersion: don't lower standards, change route or flag.)
```

## Phase 7: Common Failure Modes & Fixes

| Symptom | Diagnosis | Fix |
|---------|-----------|-----|
| TD-DFT CT state energy too low | Missing long-range exchange | Switch to CAM-B3LYP or ωB97X-D |
| TD-DFT state ordering changes with functional | Multi-reference character likely | Run CASSCF to confirm |
| CASSCF doesn't converge | Bad initial orbitals or active space | Use DFT NOs; try AVAS; adjust active space |
| CASSCF excited state higher than expected | Missing dynamic correlation | Add NEVPT2 |
| NEVPT2 correction > 1.0 eV per state | Active space too small | Enlarge active space; check NOONs |
| CASSCF NOONs all near 0 or 2 | Active space too large (dead orbitals) | Remove dead orbitals, shrink CAS |
| SA-CASSCF root flipping | States swap during optimization | Increase n_states; check orbital character |
| Negative excitation energy | Wrong reference state or broken calculation | Check ground-state stability; re-examine active space |
| Memory error during NEVPT2 | Basis too large for CAS size | Reduce basis for CAS; use RI if available |

## Checklist Before Reporting Excited-State Results

- [ ] **Ground state converged**: SCF passed all qchem-dft gates before excited-state calc
- [ ] **Method fully specified**: TD-DFT functional or CASSCF(n,m)/NEVPT2 + basis + software + version
- [ ] **Excitation type stated**: vertical absorption, vertical emission, adiabatic, 0-0
- [ ] **For TD-DFT**: NTO analysis done; CT diagnostic checked; functional appropriate for excitation type
- [ ] **For CASSCF**: Active space justified; NOONs reported; state-averaging described
- [ ] **For NEVPT2**: PT2 correction magnitude reasonable (< 1.0 eV typical)
- [ ] **Compared to reference**: Experimental λ_max / emission / literature CASPT2 values
- [ ] **Solvent effects considered**: If comparing to solution experiment, PCM included
- [ ] **Resource assessment done**: Method feasible for system size and available hardware
- [ ] **Limitations stated**: Known deficiencies of chosen method for this system type

## Integration with Vault

- **Depends on**: `qchem-dft` (ground-state geometry and SCF must pass all gates first)
- **Literature**: Use `chem-literature` for benchmark data and method comparison papers
- **Cross-agent**: CE provides candidate molecules → QE computes excited-state properties (optical gap, absorption λ) → CE uses for screening (e.g., OPV material selection)
- **Save calculations** to `research/qchem/excited-state/<date>-<system>/`
- **Git commit**: `cd <repo_root> && git add -A && git commit -m "qchem-excited-state: <description>"`

