#!/usr/bin/env python
"""Frequency + ZPE calculation from existing optimized geometries.

Reads aspirin_b97-d_def2-{svp,tzvp}_opt.xyz, runs Hessian, writes summary.json.
Skips geometry optimization (already done).
"""
import json
from pathlib import Path
from typing import List, Dict, Any

import numpy as np
from pyscf import __version__ as pyscf_version
from pyscf import gto, dft
from pyscf.hessian.thermo import harmonic_analysis

WORKDIR = Path(__file__).resolve().parent
NAME = "aspirin"
SMILES = "CC(=O)Oc1ccccc1C(=O)O"
FUNCTIONAL = "b97-d"
DISPERSION_LABEL = "built-in dispersion (XC: b97-d)"
GRID_LEVEL = 4
SCF_CONV_TOL = 1e-9
CHARGE = 0
SPIN = 0

HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KCAL_MOL = 627.509474
h = 6.62607015e-34
c = 2.99792458e10
Na = 6.02214076e23
Eh_J = 4.3597447222071e-18


def xyz_to_atom_str(xyz_text: str) -> str:
    lines = [ln.strip() for ln in xyz_text.strip().splitlines()]
    body = lines[2:]
    atom_lines = []
    for ln in body:
        if not ln:
            continue
        sym, x, y, z = ln.split()[:4]
        atom_lines.append(f"{sym} {x} {y} {z}")
    return "; ".join(atom_lines)


def zpe_from_freqs(freqs_cm1: List[float]) -> float:
    real = [f for f in freqs_cm1 if f > 0]
    zpe_Jmol = 0.5 * h * c * sum(real) * Na
    zpe_Ha = zpe_Jmol / (Eh_J * Na)
    return float(zpe_Ha)


def run_freq(xyz_path: Path, basis: str) -> Dict[str, Any]:
    print(f"\n{'='*60}")
    print(f"  Frequency analysis: {FUNCTIONAL}/{basis}")
    print(f"  Geometry from: {xyz_path.name}")
    print(f"{'='*60}\n")

    atom_str = xyz_to_atom_str(xyz_path.read_text())
    mol = gto.M(atom=atom_str, unit="Angstrom", basis=basis,
                charge=CHARGE, spin=SPIN, verbose=4)

    mf = dft.RKS(mol)
    mf.xc = FUNCTIONAL
    mf.grids.level = GRID_LEVEL
    mf.conv_tol = SCF_CONV_TOL
    mf.max_cycle = 200
    e = mf.kernel()

    if not mf.converged:
        raise RuntimeError(f"SCF not converged for {basis} (gate fail)")
    print(f"\n✅ SCF converged: E = {e:.10f} Ha")

    # Properties
    mo_e = np.array(mf.mo_energy)
    mo_occ = np.array(mf.mo_occ)
    occ_idx = np.where(mo_occ > 0)[0]
    vir_idx = np.where(mo_occ == 0)[0]
    homo = float(mo_e[occ_idx.max()])
    lumo = float(mo_e[vir_idx.min()])
    gap = lumo - homo
    dip = mf.dip_moment(unit="Debye")
    dip_norm = float(np.linalg.norm(dip))

    # Hessian + frequencies
    print(f"\nStarting Hessian for {basis} (this will take a while)...")
    hess = mf.Hessian().kernel()
    ha = harmonic_analysis(mol, hess)
    freqs = [float(x) for x in ha["freq_wavenumber"]]
    n_imag = int(sum(1 for f in freqs if f < 0))
    min_freq = float(min(freqs))

    zpe_ha = zpe_from_freqs(freqs)

    print(f"\n✅ Frequency analysis complete for {basis}")
    print(f"   n_imag = {n_imag}, min_freq = {min_freq:.2f} cm⁻¹")
    print(f"   ZPE = {zpe_ha:.8f} Ha = {zpe_ha * HARTREE_TO_KCAL_MOL:.2f} kcal/mol")

    if n_imag != 0:
        print(f"   ⚠️ WARNING: n_imag={n_imag}, frequency gate FAILED")

    return {
        "functional": FUNCTIONAL,
        "basis": basis,
        "dispersion": DISPERSION_LABEL,
        "pyscf_version": pyscf_version,
        "scf_converged": True,
        "e_tot_ha": float(e),
        "homo_ha": homo,
        "lumo_ha": lumo,
        "gap_ha": float(gap),
        "gap_ev": float(gap * HARTREE_TO_EV),
        "dipole_debye": dip_norm,
        "n_imag": n_imag,
        "min_freq_cm1": min_freq,
        "zpe_ha": zpe_ha,
        "zpe_kcal_mol": zpe_ha * HARTREE_TO_KCAL_MOL,
        "freqs_cm1": freqs,
        "geom_xyz_path": str(xyz_path),
    }


def main():
    svp_xyz = WORKDIR / "aspirin_b97-d_def2-svp_opt.xyz"
    tzvp_xyz = WORKDIR / "aspirin_b97-d_def2-tzvp_opt.xyz"

    results = []

    # --- SVP ---
    print("=" * 60)
    print("  STEP 1/2: def2-SVP frequency + ZPE")
    print("=" * 60)
    r_svp = run_freq(svp_xyz, "def2-svp")
    results.append(r_svp)

    # Write intermediate summary (in case TZVP gets interrupted)
    intermediate = {
        "molecule": NAME, "smiles": SMILES,
        "status": "svp_complete_tzvp_pending",
        "results": results,
    }
    (WORKDIR / "summary.json").write_text(json.dumps(intermediate, indent=2))
    print("\n📝 Intermediate summary.json written (SVP done)\n")

    # --- TZVP ---
    print("=" * 60)
    print("  STEP 2/2: def2-TZVP frequency + ZPE")
    print("=" * 60)
    r_tzvp = run_freq(tzvp_xyz, "def2-tzvp")
    results.append(r_tzvp)

    # Final summary
    summary = {
        "molecule": NAME,
        "smiles": SMILES,
        "charge": CHARGE,
        "spin": SPIN,
        "functional": FUNCTIONAL,
        "dispersion": DISPERSION_LABEL,
        "software": {"pyscf": pyscf_version},
        "grid_level": GRID_LEVEL,
        "scf_conv_tol": SCF_CONV_TOL,
        "status": "complete",
        "results": results,
        "basis_comparison": {
            "gap_ev(def2-svp)": r_svp["gap_ev"],
            "gap_ev(def2-tzvp)": r_tzvp["gap_ev"],
            "dipole_debye(def2-svp)": r_svp["dipole_debye"],
            "dipole_debye(def2-tzvp)": r_tzvp["dipole_debye"],
            "zpe_kcal_mol(def2-svp)": r_svp["zpe_kcal_mol"],
            "zpe_kcal_mol(def2-tzvp)": r_tzvp["zpe_kcal_mol"],
        },
    }
    (WORKDIR / "summary.json").write_text(json.dumps(summary, indent=2))

    print("\n=== Aspirin B97-D DFT summary ===\n")
    for r in results:
        tag = f"{r['functional']}/{r['basis']} ({r['dispersion']}), PySCF {r['pyscf_version']}"
        print(f"[{tag}]")
        print(f"  E_total: {r['e_tot_ha']:.10f} Ha")
        print(f"  HOMO:    {r['homo_ha']:.6f} Ha")
        print(f"  LUMO:    {r['lumo_ha']:.6f} Ha")
        print(f"  Gap:     {r['gap_ha']:.6f} Ha = {r['gap_ev']:.3f} eV")
        print(f"  Dipole:  {r['dipole_debye']:.3f} Debye")
        print(f"  ZPE:     {r['zpe_ha']:.8f} Ha = {r['zpe_kcal_mol']:.2f} kcal/mol")
        print(f"  Freq gate: n_imag={r['n_imag']} (min freq {r['min_freq_cm1']:.2f} cm⁻¹)")
        print()

    print(f"✅ All done. Summary written to: {WORKDIR / 'summary.json'}")


if __name__ == "__main__":
    main()
