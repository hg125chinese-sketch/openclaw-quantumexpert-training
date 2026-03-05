#!/usr/bin/env python
"""End-to-end DFT workflow for aspirin (acetylsalicylic acid) using PySCF.

Workflow (matches qchem-dft gates):
- SMILES -> RDKit 3D
- DFT geometry optimization with geomeTRIC (must converge)
- Single-point properties at optimized geometry: E, HOMO/LUMO/gap, dipole
- Frequency analysis at optimized geometry (gate: n_imag == 0)
- ZPE from harmonic frequencies
- Basis comparison: def2-SVP -> def2-TZVP (report properties; avoid interpreting absolute energies as convergence)

Method label (this run):
- XC: b97-d (dispersion-inclusive functional via libxc)
- Bases: def2-SVP, def2-TZVP
- Software: PySCF

Outputs are written under WORKDIR.
"""

from __future__ import annotations

import json
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
WORKDIR = Path("research/aspirin_b97d_dft")

CHARGE = 0
SPIN = 0  # 2S

FUNCTIONAL = "b97-d"  # dispersion-inclusive; supported in PySCF/libxc
DISPERSION_LABEL = "built-in dispersion (XC: b97-d)"
GRID_LEVEL = 4
SCF_CONV_TOL = 1e-9

HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KCAL_MOL = 627.509474

# constants for ZPE
h = 6.62607015e-34
c = 2.99792458e10  # cm/s
Na = 6.02214076e23
Eh_J = 4.3597447222071e-18  # J


@dataclass
class RunResult:
    functional: str
    basis: str
    dispersion: str
    pyscf_version: str

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

    freqs_cm1: List[float]
    geom_xyz_path: str


def smiles_to_xyz(smiles: str, seed: int = 7) -> str:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    ok = AllChem.EmbedMolecule(mol, params)
    if ok != 0:
        raise RuntimeError(f"RDKit EmbedMolecule failed (code {ok})")
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)

    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), f"{NAME} from SMILES {smiles}"]
    for a in mol.GetAtoms():
        pos = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():2s} {pos.x: .8f} {pos.y: .8f} {pos.z: .8f}")
    return "\n".join(lines) + "\n"


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


def build_mf(atom_str: str, basis: str) -> Tuple[gto.Mole, dft.rks.RKS]:
    mol = gto.M(atom=atom_str, unit="Angstrom", basis=basis, charge=CHARGE, spin=SPIN, verbose=4)
    mf = dft.RKS(mol)
    mf.xc = FUNCTIONAL
    mf.grids.level = GRID_LEVEL
    mf.conv_tol = SCF_CONV_TOL
    mf.max_cycle = 200
    return mol, mf


def scf_props(mf: dft.rks.RKS) -> Dict[str, Any]:
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError("SCF did not converge (gate fail)")

    mo_e = np.array(mf.mo_energy)
    mo_occ = np.array(mf.mo_occ)
    occ_idx = np.where(mo_occ > 0)[0]
    vir_idx = np.where(mo_occ == 0)[0]
    homo = float(mo_e[occ_idx.max()])
    lumo = float(mo_e[vir_idx.min()])
    gap = lumo - homo

    dip = mf.dip_moment(unit="Debye")
    dip_norm = float(np.linalg.norm(dip))

    return {
        "e_tot_ha": float(e),
        "homo_ha": homo,
        "lumo_ha": lumo,
        "gap_ha": float(gap),
        "gap_ev": float(gap * HARTREE_TO_EV),
        "dipole_debye": dip_norm,
    }


def zpe_from_freqs(freqs_cm1: List[float]) -> float:
    real = [f for f in freqs_cm1 if f > 0]
    zpe_Jmol = 0.5 * h * c * sum(real) * Na
    zpe_Ha = zpe_Jmol / (Eh_J * Na)
    return float(zpe_Ha)


def freq_and_zpe(atom_str: str, basis: str) -> Dict[str, Any]:
    mol, mf = build_mf(atom_str, basis)
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("SCF not converged before Hessian (gate fail)")

    hess = mf.Hessian().kernel()
    ha = harmonic_analysis(mol, hess)
    freqs = [float(x) for x in ha["freq_wavenumber"]]
    n_imag = int(sum(1 for f in freqs if f < 0))
    min_freq = float(min(freqs))

    zpe_ha = zpe_from_freqs(freqs)
    return {
        "freqs_cm1": freqs,
        "n_imag": n_imag,
        "min_freq_cm1": min_freq,
        "zpe_ha": zpe_ha,
        "zpe_kcal_mol": zpe_ha * HARTREE_TO_KCAL_MOL,
    }


def write_xyz(path: Path, mol: gto.Mole, comment: str) -> None:
    coords = mol.atom_coords()  # Bohr
    # convert Bohr -> Angstrom
    bohr_to_ang = 0.529177210903
    coords_ang = coords * bohr_to_ang
    lines = [str(mol.natm), comment]
    for i in range(mol.natm):
        sym = mol.atom_symbol(i)
        x, y, z = coords_ang[i]
        lines.append(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}")
    path.write_text("\n".join(lines) + "\n")


def optimize_geometry(atom_str: str, basis: str, xyz_path: Path) -> str:
    _, mf = build_mf(atom_str, basis)
    mol_opt = optimize(mf, maxsteps=200)  # gate: must converge

    # gate: SCF converged at final geometry
    mf_final = dft.RKS(mol_opt)
    mf_final.xc = FUNCTIONAL
    mf_final.grids.level = GRID_LEVEL
    mf_final.conv_tol = SCF_CONV_TOL
    mf_final.kernel()
    if not mf_final.converged:
        raise RuntimeError("SCF not converged at optimized geometry (gate fail)")

    write_xyz(
        xyz_path,
        mol_opt,
        f"Optimized {NAME} at {FUNCTIONAL}/{basis} ({DISPERSION_LABEL}), PySCF {pyscf_version}",
    )
    return xyz_path.read_text()


def run_one_basis(atom_str_init: str, basis: str) -> RunResult:
    xyz_path = WORKDIR / f"{NAME}_{FUNCTIONAL}_{basis}_opt.xyz"
    optimize_geometry(atom_str_init, basis, xyz_path)

    atom_str_opt = xyz_to_atom_str(xyz_path.read_text())
    _, mf = build_mf(atom_str_opt, basis)
    props = scf_props(mf)

    freq = freq_and_zpe(atom_str_opt, basis)
    if freq["n_imag"] != 0:
        raise RuntimeError(f"Frequency gate fail: n_imag={freq['n_imag']} (min freq {freq['min_freq_cm1']:.2f} cm^-1)")

    return RunResult(
        functional=FUNCTIONAL,
        basis=basis,
        dispersion=DISPERSION_LABEL,
        pyscf_version=pyscf_version,
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
        freqs_cm1=freq["freqs_cm1"],
        geom_xyz_path=str(xyz_path),
    )


def main() -> None:
    WORKDIR.mkdir(parents=True, exist_ok=True)

    init_xyz = smiles_to_xyz(SMILES)
    init_xyz_path = WORKDIR / f"{NAME}_rdkit_init.xyz"
    init_xyz_path.write_text(init_xyz)
    atom_init = xyz_to_atom_str(init_xyz)

    # SVP first
    r_svp = run_one_basis(atom_init, "def2-svp")

    # TZVP start from SVP optimized geometry
    atom_from_svp = xyz_to_atom_str(Path(r_svp.geom_xyz_path).read_text())
    r_tzvp = run_one_basis(atom_from_svp, "def2-tzvp")

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
        "results": [asdict(r_svp), asdict(r_tzvp)],
        "basis_comparison": {
            "gap_ev(def2-svp)": r_svp.gap_ev,
            "gap_ev(def2-tzvp)": r_tzvp.gap_ev,
            "dipole_debye(def2-svp)": r_svp.dipole_debye,
            "dipole_debye(def2-tzvp)": r_tzvp.dipole_debye,
            "zpe_kcal_mol(def2-svp)": r_svp.zpe_kcal_mol,
            "zpe_kcal_mol(def2-tzvp)": r_tzvp.zpe_kcal_mol,
        },
    }

    (WORKDIR / "summary.json").write_text(json.dumps(summary, indent=2))

    print("\n=== Aspirin DFT summary ===\n")
    for r in [r_svp, r_tzvp]:
        tag = f"{r.functional}/{r.basis} ({r.dispersion}), PySCF {r.pyscf_version}"
        print(f"[{tag}]")
        print(f"  E_total: {r.e_tot_ha:.10f} Ha")
        print(f"  HOMO:    {r.homo_ha:.6f} Ha")
        print(f"  LUMO:    {r.lumo_ha:.6f} Ha")
        print(f"  Gap:     {r.gap_ha:.6f} Ha = {r.gap_ev:.3f} eV")
        print(f"  Dipole:  {r.dipole_debye:.3f} Debye")
        print(f"  ZPE:     {r.zpe_ha:.8f} Ha = {r.zpe_kcal_mol:.2f} kcal/mol")
        print(f"  Freq gate: n_imag={r.n_imag} (min freq {r.min_freq_cm1:.2f} cm^-1)")
        print()

    print(f"Wrote outputs under: {WORKDIR}")


if __name__ == "__main__":
    main()
