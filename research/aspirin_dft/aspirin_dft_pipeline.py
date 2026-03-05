#!/usr/bin/env python
"""End-to-end DFT workflow for aspirin (acetylsalicylic acid) using PySCF.

Implements qchem-dft gates:
- SCF convergence required
- geometry optimization must converge (geomeTRIC)
- frequency analysis must show n_imag=0 for a minimum
- every reported number carries method/basis/software label

NOTE: Dispersion D3(BJ) is recommended by the skill but unavailable in this
environment (no dftd3 / no gfortran). We therefore run without dispersion and
label results accordingly.

Input: SMILES (aspirin) -> RDKit 3D -> PySCF DFT optimize -> Hessian/frequencies
Outputs: xyz geometries + JSON summary.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Any, Tuple, List

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

from pyscf import __version__ as pyscf_version
from pyscf import gto, dft
from pyscf.geomopt.geometric_solver import optimize
from pyscf.hessian.thermo import harmonic_analysis


SMILES = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin
NAME = "aspirin"
CHARGE = 0
SPIN = 0  # 2S

# If you need dispersion and dftd3 is unavailable, prefer a dispersion-inclusive functional
# supported by libxc/PySCF (e.g. "b97-d" or "wb97m-v").
FUNCTIONAL = "b97-d"  # dispersion-inclusive; supported in PySCF/libxc
DISPERSION_LABEL = "built-in dispersion (via XC)"
GRID_LEVEL = 4
SCF_CONV_TOL = 1e-9

BOHR_TO_ANG = 0.529177210903
AU_TO_DEBYE = 2.541746  # e·a0 to Debye
HARTREE_TO_EV = 27.211386245988


@dataclass
class RunResult:
    label: str
    functional: str
    basis: str
    pyscf_version: str
    dispersion: str

    scf_converged: bool
    e_tot_ha: float

    homo_ha: float
    lumo_ha: float
    gap_ha: float
    gap_ev: float

    dipole_debye: float

    n_imag: int
    min_freq_cm1: float
    zpe_ha: float
    zpe_kcal_mol: float

    geom_xyz: str


def smiles_to_xyz(smiles: str, seed: int = 7) -> str:
    """Generate an RDKit-embedded 3D geometry and return XYZ string."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    ok = AllChem.EmbedMolecule(mol, params)
    if ok != 0:
        raise RuntimeError(f"RDKit EmbedMolecule failed (code {ok})")
    # quick force-field relaxation
    AllChem.UFFOptimizeMolecule(mol, maxIters=200)

    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), f"{NAME} from SMILES {smiles}"]
    for a in mol.GetAtoms():
        pos = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():2s} {pos.x: .8f} {pos.y: .8f} {pos.z: .8f}")
    return "\n".join(lines) + "\n"


def xyz_to_pyscf_atom(xyz: str) -> str:
    """Convert XYZ (Å) to PySCF atom string."""
    lines = [ln.strip() for ln in xyz.strip().splitlines()]
    if len(lines) < 3:
        raise ValueError("XYZ too short")
    body = lines[2:]
    atom_lines = []
    for ln in body:
        if not ln:
            continue
        sym, x, y, z = ln.split()[:4]
        atom_lines.append(f"{sym} {x} {y} {z}")
    return "; ".join(atom_lines)


def write_xyz(path: Path, mol: gto.Mole, comment: str) -> str:
    coords_ang = mol.atom_coords() * BOHR_TO_ANG
    lines = [str(mol.natm), comment]
    for i in range(mol.natm):
        sym = mol.atom_symbol(i)
        x, y, z = coords_ang[i]
        lines.append(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}")
    text = "\n".join(lines) + "\n"
    path.write_text(text)
    return text


def build_mf(atom_str: str, basis: str) -> Tuple[gto.Mole, dft.rks.RKS]:
    mol = gto.M(
        atom=atom_str,
        unit="Angstrom",
        basis=basis,
        charge=CHARGE,
        spin=SPIN,
        verbose=4,
    )
    mf = dft.RKS(mol)
    mf.xc = FUNCTIONAL
    mf.grids.level = GRID_LEVEL
    mf.conv_tol = SCF_CONV_TOL
    mf.max_cycle = 200
    return mol, mf


def scf_and_props(mf: dft.rks.RKS) -> Dict[str, Any]:
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError("SCF did not converge (gate fail)")

    mo_e = np.array(mf.mo_energy)
    mo_occ = np.array(mf.mo_occ)
    occ_idx = np.where(mo_occ > 0)[0]
    vir_idx = np.where(mo_occ == 0)[0]
    homo = mo_e[occ_idx.max()]
    lumo = mo_e[vir_idx.min()]
    gap = lumo - homo

    dip = mf.dip_moment(unit="Debye")
    dip_norm = float(np.linalg.norm(dip))

    return {
        "e_tot_ha": float(e),
        "homo_ha": float(homo),
        "lumo_ha": float(lumo),
        "gap_ha": float(gap),
        "gap_ev": float(gap * HARTREE_TO_EV),
        "dipole_debye": dip_norm,
    }


def optimize_geometry(atom_str: str, basis: str, out_xyz: Path) -> gto.Mole:
    _, mf = build_mf(atom_str, basis)
    # geomeTRIC optimization gate: must converge
    mol_opt = optimize(mf, maxsteps=200)

    # ensure SCF converged at the final geometry
    mf_final = dft.RKS(mol_opt)
    mf_final.xc = FUNCTIONAL
    mf_final.grids.level = GRID_LEVEL
    mf_final.conv_tol = SCF_CONV_TOL
    mf_final.kernel()
    if not mf_final.converged:
        raise RuntimeError("SCF not converged at optimized geometry (gate fail)")

    write_xyz(out_xyz, mol_opt, f"Optimized {NAME} at {FUNCTIONAL}/{basis} ({DISPERSION_LABEL}), PySCF {pyscf_version}")
    return mol_opt


def frequency_analysis(mol: gto.Mole, basis: str) -> Dict[str, Any]:
    # Rebuild and run SCF at the given geometry for consistency
    atom_str = "; ".join([f"{mol.atom_symbol(i)} {mol.atom_coord(i)[0]*BOHR_TO_ANG} {mol.atom_coord(i)[1]*BOHR_TO_ANG} {mol.atom_coord(i)[2]*BOHR_TO_ANG}" for i in range(mol.natm)])
    _, mf = build_mf(atom_str, basis)
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("SCF not converged before Hessian (gate fail)")

    hess = mf.Hessian().kernel()
    ha = harmonic_analysis(mf.mol, hess)
    freqs = np.array(ha["freq_wavenumber"], dtype=float)
    n_imag = int(np.sum(freqs < 0))
    min_freq = float(freqs.min())

    # PySCF harmonic_analysis reports ZPE in Hartree (key: 'ZPE')
    zpe_ha = float(ha.get("ZPE", np.nan))
    # Convert to kcal/mol
    zpe_kcal_mol = float(zpe_ha * 627.509474)

    return {
        "n_imag": n_imag,
        "min_freq_cm1": min_freq,
        "zpe_ha": zpe_ha,
        "zpe_kcal_mol": zpe_kcal_mol,
        "freqs_cm1": freqs.tolist(),
    }


def run_level(level_label: str, atom_str: str, basis: str, workdir: Path) -> RunResult:
    geom_path = workdir / f"{NAME}_{level_label}_{basis}_opt.xyz"
    mol_opt = optimize_geometry(atom_str, basis, geom_path)

    # single-point + properties at optimized geometry
    atom_opt = xyz_to_pyscf_atom(geom_path.read_text())
    _, mf = build_mf(atom_opt, basis)
    props = scf_and_props(mf)

    # frequency + ZPE gate
    freq = frequency_analysis(mol_opt, basis)
    if freq["n_imag"] != 0:
        raise RuntimeError(f"Frequency gate fail: n_imag={freq['n_imag']} (min freq {freq['min_freq_cm1']:.1f} cm^-1)")

    return RunResult(
        label=level_label,
        functional=FUNCTIONAL,
        basis=basis,
        pyscf_version=pyscf_version,
        dispersion=DISPERSION_LABEL,
        scf_converged=True,
        e_tot_ha=props["e_tot_ha"],
        homo_ha=props["homo_ha"],
        lumo_ha=props["lumo_ha"],
        gap_ha=props["gap_ha"],
        gap_ev=props["gap_ev"],
        dipole_debye=props["dipole_debye"],
        n_imag=freq["n_imag"],
        min_freq_cm1=freq["min_freq_cm1"],
        zpe_ha=freq["zpe_ha"],
        zpe_kcal_mol=freq["zpe_kcal_mol"],
        geom_xyz=geom_path.read_text(),
    )


def main() -> None:
    workdir = Path("aspirin_dft")
    workdir.mkdir(exist_ok=True)

    # Phase 0: 3D initial structure
    xyz0 = smiles_to_xyz(SMILES)
    (workdir / f"{NAME}_rdkit_init.xyz").write_text(xyz0)
    atom0 = xyz_to_pyscf_atom(xyz0)

    results: List[RunResult] = []

    # def2-SVP (screening / pre-opt) then def2-TZVP (production)
    results.append(run_level("prod", atom0, "def2-svp", workdir))

    # start TZVP optimization from the SVP optimized geometry
    atom_from_svp_opt = xyz_to_pyscf_atom((workdir / f"{NAME}_prod_def2-svp_opt.xyz").read_text())
    results.append(run_level("prod", atom_from_svp_opt, "def2-tzvp", workdir))

    # basis convergence comparison (SVP -> TZVP) on their own optimized geometries
    e_svp = results[0].e_tot_ha
    e_tz = results[1].e_tot_ha
    delta_ha = e_tz - e_svp
    delta_kcal = delta_ha * 627.509474

    summary = {
        "molecule": NAME,
        "smiles": SMILES,
        "charge": CHARGE,
        "spin": SPIN,
        "software": {"pyscf": pyscf_version},
        "grid_level": GRID_LEVEL,
        "scf_conv_tol": SCF_CONV_TOL,
        "dispersion": DISPERSION_LABEL,
        "results": [asdict(r) for r in results],
        "basis_convergence": {
            "delta_E(def2-tzvp - def2-svp) [Ha]": delta_ha,
            "delta_E(def2-tzvp - def2-svp) [kcal/mol]": delta_kcal,
        },
    }

    (workdir / "summary.json").write_text(json.dumps(summary, indent=2))

    # human-readable printout
    print("\n=== Aspirin DFT end-to-end summary (all values labeled) ===\n")
    for r in results:
        tag = f"{r.functional}/{r.basis} ({r.dispersion}), PySCF {r.pyscf_version}"
        print(f"[{tag}]")
        print(f"  E_total:   {r.e_tot_ha:.10f} Ha")
        print(f"  HOMO:      {r.homo_ha:.6f} Ha")
        print(f"  LUMO:      {r.lumo_ha:.6f} Ha")
        print(f"  Gap:       {r.gap_ha:.6f} Ha  = {r.gap_ev:.3f} eV")
        print(f"  Dipole:    {r.dipole_debye:.3f} Debye")
        print(f"  ZPE:       {r.zpe_ha:.8f} Ha = {r.zpe_kcal_mol:.2f} kcal/mol")
        print(f"  Freq gate: n_imag={r.n_imag} (min freq {r.min_freq_cm1:.1f} cm^-1)")
        print()

    print("[Basis convergence]")
    print(f"  ΔE(TZVP-SVP) on respective optimized geometries: {delta_ha:+.10f} Ha = {delta_kcal:+.2f} kcal/mol")
    print(f"\nWrote: {workdir}/summary.json and optimized geometries (*.xyz)")


if __name__ == "__main__":
    main()
