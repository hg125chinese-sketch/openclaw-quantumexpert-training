---
name: qchem-dft
description: DFT calculation setup, functional/basis set selection, SCF convergence diagnostics, geometry optimization, frequency analysis, and result quality control. The foundational quantum chemistry skill — every other qchem skill builds on this.
homepage: https://pyscf.org
metadata: { "openclaw": { "emoji": "⚛️", "requires": { "bins": ["python3"], "python": ["pyscf", "numpy", "scipy", "matplotlib"] } } }
---

# DFT Calculations: Setup, Execution, Diagnosis

The core quantum chemistry skill. Covers functional and basis set selection, SCF convergence, geometry optimization, frequency analysis, and systematic quality control. Everything downstream (excited states, thermochemistry, benchmarking) depends on getting DFT right.

## When to Use

- User asks to compute electronic energy, geometry, or properties of a molecule
- User needs to choose a DFT functional or basis set
- User encounters SCF convergence problems
- User wants to optimize a molecular geometry
- User needs frequency analysis (IR/Raman, ZPE, thermochemistry, transition state validation)
- User asks about basis set superposition error or basis set convergence
- User wants to compare DFT methods or validate against experiment/higher-level theory
- Another qchem skill needs a geometry or single-point energy as input

## Core Philosophy

1. **SCF convergence is not optional.** If it didn't converge, nothing downstream is real. Diagnose first, always.
2. **Every number needs a method label.** Never report an energy, gap, or barrier without: functional, basis set, dispersion correction, software, version. Naked numbers are meaningless numbers.
3. **Know your approximations.** DFT is not exact. No functional is universal. Be honest about what your level of theory can and cannot resolve.
4. **Compare to reference, always.** Experimental values when available. CCSD(T)/CBS as computational gold standard. Never present results in a vacuum.
5. **Plan before computing.** Lay out: method → basis → geometry → property → QC checks. Get approval, then execute.
6. **Converge the basis before trusting the answer.** A result that changes by 5 kcal/mol when you go from DZ to TZ is not a result — it's a basis set artifact.

## Phase 1: Functional Selection

### 1.1 Decision Tree

```
What are you calculating?
│
├── Organic small molecule (main-group, no metals)
│   ├── Geometry / thermochemistry → B3LYP-D3(BJ) / def2-TZVP (safe default)
│   ├── Reaction barrier → ωB97X-D or M06-2X (better for barriers)
│   ├── Non-covalent interactions → ωB97X-D or B3LYP-D3(BJ) (dispersion critical)
│   └── HOMO-LUMO gap / excitation → ωB97X-D or PBE0 (range-separated preferred for gap)
│
├── Transition metal complex
│   ├── Geometry → PBE0-D3(BJ) or TPSSh-D3(BJ) (start here)
│   ├── Spin state energetics → TPSSh / B3LYP* (HF exchange sensitivity!)
│   │   └── ALWAYS do multi-functional comparison: TPSSh / PBE0 / B3LYP / B3LYP*
│   │       (spin splitting varies by 5–20 kcal/mol across functionals)
│   └── Redox potentials → B3LYP-D3 + solvation model (CPCM/SMD)
│
├── Extended / periodic system
│   └── PBE-D3 or r²SCAN (revPBE for specific benchmarks)
│
└── Unsure / general purpose
    └── ωB97X-D / def2-TZVP (good all-rounder, rarely terrible)
```

### 1.2 Dispersion Correction: Always On

```
Dispersion correction decision:
├── D3(BJ) — default for most functionals (Grimme 2011)
├── D4 — newer, slightly better for metals, not always available
├── -D (built-in, e.g., ωB97X-D) — already includes dispersion
└── No dispersion → ONLY acceptable if you explicitly justify why
    (e.g., benchmarking against non-dispersion-corrected reference)

Rule: If your functional name doesn't end in -D, add D3(BJ).

**Environment-aware rule (mandatory):**
- If D3/D4 tooling is unavailable in your environment (cannot import/build/run), **do not silently drop dispersion**.
- Choose one of:
  1) **Switch to a dispersion-inclusive functional available in your backend**, e.g. **B97-D** or **ωB97M-V** (VV10), and document the switch.
  2) **Stop and ask for explicit approval** to run without dispersion, and justify the exception.
```

### 1.3 HF Exchange Sensitivity Warning (Transition Metals)

Spin state energetics of transition metal complexes are **critically sensitive** to the fraction of HF exchange:

| Functional | %HF exchange | Spin state bias |
|-----------|-------------|----------------|
| PBE / TPSS | 0% | Strongly favors low-spin |
| TPSSh | 10% | Mild LS preference (often best compromise) |
| B3LYP* | 15% | Roos' modified B3LYP for TM spin states |
| B3LYP | 20% | Moderate HS preference |
| PBE0 | 25% | Stronger HS preference |
| M06 | 27% | Variable; not recommended for spin states |

**Rule**: For spin state splitting, NEVER trust a single functional. Run at least TPSSh + PBE0 + B3LYP and compare. If they disagree by > 5 kcal/mol, flag as "DFT-unreliable, multi-reference may be needed."

## Phase 2: Basis Set Selection

### 2.1 The def2 Family (Recommended Default)

```
Basis set ladder (Ahlrichs/Karlsruhe):
│
├── def2-SVP      (~DZ)  — geometry pre-optimization, quick screening
│   └── Use for: initial geometry, large systems, exploratory
│
├── def2-TZVP     (~TZ)  — production geometry + single point
│   └── Use for: final geometries, reliable energetics
│
├── def2-TZVPP    (~TZ+) — improved single point
│   └── Use for: energy refinement on TZVP geometry
│
└── def2-QZVPP    (~QZ)  — near CBS, expensive
    └── Use for: benchmark, basis set convergence check

ECPs: def2 family includes ECPs automatically for elements beyond Kr.
Auxiliary bases: def2/J (Coulomb fitting), def2/JK (exchange), def2-TZVP/C (correlation)
```

### 2.2 Basis Set Convergence Test Protocol

**Always run this when reporting quantitative energetics** (barriers, binding energies, spin splitting):

```python
#!/opt/conda/envs/chem/bin/python
"""Basis set convergence test for a given molecule and property."""
from pyscf import gto, dft
import numpy as np

def basis_convergence_test(atom_str, charge=0, spin=0, functional='b3lyp'):
    """Run DFT with increasing basis sets; check energy convergence.

    Returns dict of {basis: energy} for convergence analysis.
    """
    bases = ['def2-svp', 'def2-tzvp', 'def2-tzvpp', 'def2-qzvpp']
    results = {}

    for basis in bases:
        mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin,
                     verbose=0)
        mf = dft.RKS(mol) if spin == 0 else dft.UKS(mol)
        mf.xc = functional
        mf.grids.level = 4  # fine grid
        mf.conv_tol = 1e-10

        try:
            e = mf.kernel()
            if mf.converged:
                results[basis] = e
                print(f"  {basis:15s}  E = {e:.10f} Ha  ✓ converged")
            else:
                results[basis] = None
                print(f"  {basis:15s}  ✗ NOT converged")
        except Exception as ex:
            results[basis] = None
            print(f"  {basis:15s}  ✗ FAILED: {ex}")

    # Convergence check
    energies = [v for v in results.values() if v is not None]
    if len(energies) >= 2:
        deltas = [abs(energies[i+1] - energies[i]) * 627.509
                  for i in range(len(energies)-1)]  # convert to kcal/mol
        print(f"\n  Δ(SVP→TZ):  {deltas[0]:.2f} kcal/mol" if len(deltas) > 0 else "")
        print(f"  Δ(TZ→TZP):  {deltas[1]:.2f} kcal/mol" if len(deltas) > 1 else "")
        print(f"  Δ(TZP→QZ):  {deltas[2]:.2f} kcal/mol" if len(deltas) > 2 else "")

        if deltas[-1] > 1.0:
            print("  ⚠️ WARNING: Not converged at largest basis. QZ+ or extrapolation needed.")
        elif deltas[-1] > 0.3:
            print("  ⚠️ CAUTION: Marginally converged. Report with error bar.")
        else:
            print("  ✅ Basis set converged (Δ < 0.3 kcal/mol).")

    return results
```

### 2.3 Basis Set Convergence Thresholds

| Property type | Δ(TZ→QZ) acceptable | Action if exceeded |
|--------------|---------------------|-------------------|
| Relative energy / barrier | < 0.5 kcal/mol | Use QZ or extrapolate to CBS |
| Spin state splitting | < 1.0 kcal/mol | Use QZ; if still > 1.0, flag as basis-sensitive |
| Binding energy | < 0.3 kcal/mol | Counterpoise correction + QZ |
| Geometry (bond length) | < 0.005 Å | TZ usually sufficient |
| Vibrational frequency | < 10 cm⁻¹ | TZ usually sufficient |

## Phase 3: SCF Convergence

### 3.1 Default Settings (Start Here)

```python
#!/opt/conda/envs/chem/bin/python
"""Standard DFT single-point calculation with proper defaults."""
from pyscf import gto, dft

mol = gto.M(
    atom='''
    O  0.000  0.000  0.117
    H  0.000  0.757 -0.469
    H  0.000 -0.757 -0.469
    ''',
    basis='def2-tzvp',
    charge=0,
    spin=0,       # 2S (number of unpaired electrons)
    verbose=4,    # enough detail to diagnose problems
)

mf = dft.RKS(mol)     # RKS for closed-shell; UKS for open-shell
mf.xc = 'b3lyp'       # functional
mf.grids.level = 4    # integration grid fineness (3=default, 4=fine, 5=ultrafine)
mf.conv_tol = 1e-9    # energy convergence (Ha); default 1e-9 is good
mf.max_cycle = 200    # max SCF iterations; increase for difficult cases

# Dispersion (D3BJ via dftd3 interface if available)
try:
    from pyscf import dftd3
    mf = dftd3.dftd3(mf)
    print("D3(BJ) dispersion correction enabled")
except ImportError:
    print("⚠️ dftd3 not available; running without dispersion")

e_tot = mf.kernel()

if mf.converged:
    print(f"\n✅ SCF converged: E = {e_tot:.10f} Ha")
else:
    print(f"\n❌ SCF NOT converged after {mf.max_cycle} cycles")
    print("   → See Phase 3.2 for troubleshooting")
```

### 3.2 SCF Troubleshooting Decision Tree

```
SCF not converging?
│
├── Step 1: Increase max_cycle to 500
│   └── Sometimes it just needs more iterations
│
├── Step 2: Try DIIS damping
│   mf.diis_space = 12          # increase DIIS history
│   mf.damp = 0.5               # damping factor (0–1, larger = more damping)
│   └── Good for oscillating SCF
│
├── Step 3: Level shifting
│   mf.level_shift = 0.5        # shift virtual orbitals up by 0.5 Ha
│   └── Good for small HOMO-LUMO gap / near-degenerate states
│
├── Step 4: Change initial guess
│   mf.init_guess = 'atom'      # default; try 'minao' or 'huckel'
│   # Or use converged smaller-basis result:
│   mf_small = dft.RKS(gto.M(atom=atom_str, basis='def2-svp')).run()
│   mf.init_guess = 'chkfile'
│   mf.__dict__.update(dm0=mf_small.make_rdm1())
│   └── Smaller basis → converge → project to larger basis
│
├── Step 5: Finer grid
│   mf.grids.level = 5          # ultrafine grid
│   └── Helps if oscillation is grid-related (rare but real)
│
├── Step 6: Second-order SCF (Newton)
│   mf = mf.newton()            # quadratic convergence near solution
│   └── Expensive per step but powerful for stubborn cases
│
└── Step 7: Fractional occupation (smearing)
    mf = mf.apply(dft.addons.frac_occ_)
    └── For metallic / near-degenerate systems; thermal smearing helps
        ⚠️ Final energy must be extrapolated to T=0
```

### 3.3 Open-Shell Systems: UKS vs ROKS

```
Open-shell molecule?
│
├── Doublet / high-spin state → UKS (unrestricted)
│   mf = dft.UKS(mol)
│   # ALWAYS check spin contamination:
│   print(f"<S²> = {mf.spin_square()}")
│   # Expected: S(S+1) for spin S
│   # |<S²> - S(S+1)| > 0.1 → spin contamination, results suspect
│
├── Singlet diradical / broken-symmetry → UKS with BS initial guess
│   # Generate broken-symmetry guess:
│   mf = dft.UKS(mol)
│   mf.init_guess = 'atom'    # or manually break alpha/beta symmetry
│   # Run stability analysis after convergence:
│   mf.stability()
│   └── If instability found → re-optimize with broken-symmetry solution
│
└── Want to avoid spin contamination → ROKS (restricted open-shell)
    mf = dft.ROKS(mol)
    └── Slower, sometimes harder to converge, but no spin contamination
```

### 3.4 Spin Contamination Check

```python
def check_spin_contamination(mf, threshold=0.1):
    """Check <S²> deviation from ideal value.

    For UKS/UHF, spin contamination is a common artifact.
    Large deviation means the wavefunction is mixing in higher spin states.
    """
    s2, mult = mf.spin_square()
    spin = mf.mol.spin / 2.0  # S quantum number
    s2_exact = spin * (spin + 1)
    deviation = abs(s2 - s2_exact)

    print(f"  <S²> = {s2:.4f} (exact: {s2_exact:.4f}, deviation: {deviation:.4f})")
    print(f"  2S+1 = {mult:.4f} (expected: {2*spin+1:.1f})")

    if deviation > threshold:
        print(f"  ⚠️ Spin contamination detected (Δ<S²> = {deviation:.4f} > {threshold})")
        print(f"     Possible actions:")
        print(f"     1) Switch to ROKS (restricted open-shell)")
        print(f"     2) Use broken-symmetry + spin projection")
        print(f"     3) Consider multi-reference method (CASSCF)")
        return False
    else:
        print(f"  ✅ Spin contamination acceptable")
        return True
```

## Phase 4: Geometry Optimization

### 4.1 Standard Geometry Optimization

```python
#!/opt/conda/envs/chem/bin/python
"""Geometry optimization with PySCF + geomeTRIC."""
from pyscf import gto, dft
from pyscf.geomopt.geometric_solver import optimize

mol = gto.M(
    atom='''...''',  # initial geometry (Å)
    basis='def2-tzvp',
    charge=0,
    spin=0,
)

mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.grids.level = 4

# geomeTRIC optimizer (recommended; install: pip install geometric)
mol_opt = optimize(mf, maxsteps=200)

print(f"\nOptimized geometry:")
print(mol_opt.atom_coords() * 0.529177)  # Bohr → Å

# IMPORTANT: save optimized geometry for downstream use
with open('optimized_geometry.xyz', 'w') as f:
    f.write(f"{mol_opt.natm}\n")
    f.write(f"Optimized at {mf.xc}/def2-tzvp\n")
    for i in range(mol_opt.natm):
        sym = mol_opt.atom_symbol(i)
        x, y, z = mol_opt.atom_coords()[i] * 0.529177  # Bohr → Å
        f.write(f"{sym}  {x:.6f}  {y:.6f}  {z:.6f}\n")
```

### 4.2 Optimization QC Checklist

After every geometry optimization, verify:

```
Optimization complete?
│
├── Check 1: Did it converge?
│   └── geomeTRIC reports "Converged!" → ✅
│   └── "Not converged" or max steps reached → ❌ increase maxsteps or check geometry
│
├── Check 2: Is SCF converged at final geometry?
│   └── Check mf.converged == True at the final point
│
├── Check 3: Frequency analysis (Phase 5)
│   └── All frequencies real → ✅ true minimum
│   └── One imaginary frequency → transition state (TS)
│   └── Multiple imaginary → higher-order saddle point or bad geometry
│
├── Check 4: Reasonable geometry?
│   └── Bond lengths within ±0.05 Å of experimental/reference
│   └── No unreasonably short contacts (< 1.0 Å for non-H)
│   └── Coordination number makes chemical sense
│
└── Check 5: For open-shell → spin contamination at final geometry
```

## Phase 5: Frequency Analysis

### 5.1 Harmonic Frequency Calculation

```python
#!/opt/conda/envs/chem/bin/python
"""Frequency analysis at optimized geometry."""
from pyscf import gto, dft
from pyscf.hessian import rks as rks_hess  # or uks for open-shell
import numpy as np

# Load optimized geometry
mol = gto.M(
    atom='''...''',  # use optimized geometry from Phase 4
    basis='def2-tzvp',
    charge=0,
    spin=0,
)

mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.grids.level = 4
mf.kernel()

# Hessian → frequencies
hess = mf.Hessian().kernel()

# Convert Hessian to frequencies (PySCF provides utilities)
from pyscf.hessian.thermo import harmonic_analysis
results = harmonic_analysis(mol, hess)
freqs = results['freq_wavenumber']  # in cm⁻¹

print("\nVibrational frequencies (cm⁻¹):")
for i, f in enumerate(freqs):
    marker = "⚠️ IMAGINARY" if f < 0 else ""
    print(f"  Mode {i+1:3d}: {f:8.1f}  {marker}")

n_imag = sum(1 for f in freqs if f < 0)
print(f"\nImaginary frequencies: {n_imag}")
if n_imag == 0:
    print("✅ True minimum confirmed")
elif n_imag == 1:
    print("⚠️ Transition state (1 imaginary frequency)")
    print("   If you wanted a minimum, see Phase 5.2")
else:
    print("❌ Higher-order saddle point — geometry is problematic")
```

### 5.2 Imaginary Frequency Troubleshooting

```
Imaginary frequency found at "optimized" geometry?
│
├── Small imaginary (> -50 cm⁻¹, e.g., -30 to -50)?
│   ├── Often numerical noise from:
│   │   ├── Coarse integration grid → increase grids.level to 5
│   │   ├── Loose optimization convergence → tighten thresholds
│   │   └── Floppy torsion / methyl rotation → may be unavoidable
│   ├── Fix attempts:
│   │   ├── 1. Re-optimize with tighter convergence + finer grid
│   │   ├── 2. Displace along imaginary mode (±0.05 Å) → re-optimize
│   │   └── 3. If floppy mode persists → document it, note effect on thermo
│   └── Report: "Frequency X cm⁻¹ identified as [torsion/rotation],
│        treated as [re-optimized / documented soft mode]"
│
├── Medium imaginary (-50 to -200 cm⁻¹)?
│   └── Likely a genuine saddle point
│   ├── Visualize the mode (what atoms move, what coordinate?)
│   ├── Displace along mode in BOTH directions → re-optimize both
│   ├── Take the lower-energy structure → re-check frequencies
│   └── May indicate a different conformer is the true minimum
│
└── Large imaginary (< -200 cm⁻¹)?
    └── Definitely a saddle point or wrong electronic state
    ├── Check: is the electronic state correct? (spin, charge, symmetry)
    ├── Check: did optimization actually converge? (sometimes false convergence)
    └── Start from a different initial geometry
```

### 5.3 Thermochemistry from Frequencies

```python
def compute_thermochemistry(mol, freqs, e_elec, T=298.15, P=101325.0):
    """Compute thermodynamic quantities from harmonic frequencies.

    Args:
        mol: PySCF Mole object (for mass and geometry)
        freqs: frequencies in cm⁻¹ (real only; imaginary excluded)
        e_elec: electronic energy in Hartree
        T: temperature in K
        P: pressure in Pa

    Returns:
        dict with E_elec, ZPE, H, G in Hartree and kcal/mol
    """
    import numpy as np

    # Constants
    h = 6.62607015e-34    # J·s
    kB = 1.380649e-23     # J/K
    c = 2.99792458e10     # cm/s
    Na = 6.02214076e23    # mol⁻¹
    Ha_to_kcal = 627.509474

    # Filter real frequencies (exclude imaginary and < 50 cm⁻¹ for RRHO issues)
    real_freqs = [f for f in freqs if f > 0]

    # ZPE
    zpe_J = 0.5 * h * c * sum(real_freqs) * Na  # J/mol
    zpe_Ha = zpe_J / (Ha_to_kcal * 4184)

    # Vibrational contributions (harmonic oscillator)
    u_vib = 0.0
    s_vib = 0.0
    for nu in real_freqs:
        x = h * c * nu / (kB * T)
        u_vib += h * c * nu * Na * (0.5 + 1.0 / (np.exp(x) - 1))  # J/mol
        s_vib += kB * Na * (x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)))  # J/(mol·K)

    H = e_elec + zpe_Ha  # simplified; full H includes translational + rotational
    G = H - T * s_vib / (Ha_to_kcal * 4184)  # simplified

    results = {
        'E_elec (Ha)': e_elec,
        'ZPE (kcal/mol)': zpe_Ha * Ha_to_kcal,
        'H (Ha)': H,
        'G (Ha)': G,
        'T (K)': T,
        'n_real_freqs': len(real_freqs),
    }

    print(f"\n=== Thermochemistry at {T:.1f} K ===")
    for k, v in results.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.6f}")
        else:
            print(f"  {k}: {v}")

    return results
```

## Phase 6: Common Failure Modes & Fixes

| Symptom | Diagnosis | Fix |
|---------|-----------|-----|
| SCF oscillates, never converges | HOMO-LUMO near-degeneracy or bad initial guess | Level shift (0.3–1.0 Ha); DIIS damping; smaller basis first |
| SCF converges to wrong state | Multiple local minima on SCF surface | Stability analysis; try different init_guess; broken-symmetry |
| Geometry optimization oscillates | Flat PES or competing minima | Tighter grid; change optimizer (geomeTRIC → BFGS); check for symmetry breaking |
| Imaginary frequency at "minimum" | Not a true minimum | See Phase 5.2 |
| Spin contamination (⟨S²⟩ wrong) | UKS mixing spin states | ROKS; broken-symmetry + projection; multi-reference |
| Energy changes > 5 kcal/mol with basis | Basis set incomplete | Go up the ladder (SVP → TZ → QZ); check Phase 2.2 |
| Dispersion-bound complex falls apart | Missing dispersion correction | Add D3(BJ) or D4; this is mandatory for non-covalent interactions |
| Transition metal spin state "wrong" | Functional HF exchange dependence | Multi-functional comparison (Phase 1.3); consider CASSCF |
| NaN / infinite energy | Basis set linear dependency or broken geometry | Check geometry; remove diffuse functions; increase linear dependency threshold |

## Phase 7: Solvation Models

### 7.1 When to Use Solvation

```
Include solvation when:
├── Comparing to solution-phase experiments (always)
├── Computing redox potentials or pKa
├── Charged species (large solvation energy)
├── Conformational preferences that depend on polarity
└── Reaction barriers in solution

Skip solvation when:
├── Gas-phase benchmark (comparing to gas-phase experiment)
├── Solid-state / periodic calculations
└── Initial exploration (add solvation at refinement stage)
```

### 7.2 PySCF Solvation Setup

```python
from pyscf import solvent

mf = dft.RKS(mol)
mf.xc = 'b3lyp'

# PCM (Polarizable Continuum Model)
mf = solvent.ddCOSMO(mf)  # ddCOSMO is PySCF's default PCM variant
mf.with_solvent.eps = 78.39  # water dielectric constant

# Common solvents:
# Water: eps=78.39  |  DMSO: eps=46.83  |  DCM: eps=8.93
# THF: eps=7.58     |  Toluene: eps=2.38 |  Hexane: eps=1.88

e_sol = mf.kernel()
```

## Checklist Before Reporting

- [ ] **Method fully specified**: Functional + dispersion + basis set + software + version
- [ ] **SCF converged**: `mf.converged == True`; convergence threshold stated
- [ ] **Spin contamination checked** (if open-shell): ⟨S²⟩ deviation < 0.1
- [ ] **Basis set convergence tested** (if quantitative energetics): Δ(TZ→QZ) within thresholds
- [ ] **Geometry verified**: Optimization converged; bond lengths reasonable
- [ ] **Frequencies checked** (if geometry optimized): No unexpected imaginary frequencies
- [ ] **Compared to reference**: Experimental value and/or higher-level calculation cited
- [ ] **Solvation included** (if comparing to solution experiment): Model stated
- [ ] **Units explicit**: Ha, eV, kcal/mol, kJ/mol — always labeled, never assumed
- [ ] **Reproducibility**: Input geometry saved; all parameters logged; script archived

## Unit Conversion Reference

```
1 Hartree  = 627.509474 kcal/mol
           = 2625.4996 kJ/mol
           = 27.211386 eV
           = 219474.63 cm⁻¹

1 eV       = 23.0605 kcal/mol
           = 96.4853 kJ/mol
           = 8065.54 cm⁻¹

1 Bohr     = 0.529177 Å
1 Å        = 1.889726 Bohr
```

## Integration with Vault

- **Literature**: Use `chem-literature` for method validation papers and experimental reference data
- **ADMET**: After QE provides electronic properties, CE can use `chem-admet` for drug-likeness
- **Docking**: QE-computed partial charges or ESP can feed into `chem-docking` for better scoring
- **Molecular generation**: CE's `chem-molgen` generates candidates → QE validates electronic structure
- **Save calculations** to `research/qchem/calculations/<date>-<system>/`
- **Git commit**: `cd /home/node/.openclaw/workspace-quantumexpert && git add -A && git commit -m "qchem-dft: <description>"`
