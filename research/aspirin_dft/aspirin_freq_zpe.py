#!/usr/bin/env python
"""Compute harmonic frequencies and ZPE for aspirin at given method/basis.

Reads optimized XYZ geometries written by aspirin_dft_pipeline.py and recomputes:
- SCF energy (gate: converged)
- Harmonic frequencies via analytic Hessian (RKS)
- ZPE from frequencies

Outputs: updates aspirin_dft/summary.json with freqs_cm1, zpe_ha, zpe_kcal_mol, min_freq_cm1, n_imag.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pyscf import gto, dft
from pyscf import __version__ as pyscf_version
from pyscf.hessian.thermo import harmonic_analysis

NAME = "aspirin"
WORKDIR = Path("aspirin_dft")
SUMMARY_PATH = WORKDIR / "summary.json"

FUNCTIONAL = "pbe0"
GRID_LEVEL = 4
SCF_CONV_TOL = 1e-9
CHARGE = 0
SPIN = 0

HARTREE_TO_KCAL_MOL = 627.509474

# constants
h = 6.62607015e-34
c = 2.99792458e10  # cm/s
Na = 6.02214076e23
Eh_J = 4.3597447222071e-18  # J


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


def zpe_from_freqs_cm1(freqs_cm1: list[float]) -> float:
    real = [f for f in freqs_cm1 if f > 0]
    zpe_Jmol = 0.5 * h * c * sum(real) * Na
    zpe_Ha = zpe_Jmol / (Eh_J * Na)
    return float(zpe_Ha)


def run_for_basis(basis: str, xyz_path: Path) -> dict:
    atom_str = xyz_to_atom_str(xyz_path.read_text())
    mol = gto.M(atom=atom_str, unit="Angstrom", basis=basis, charge=CHARGE, spin=SPIN, verbose=4)
    mf = dft.RKS(mol)
    mf.xc = FUNCTIONAL
    mf.grids.level = GRID_LEVEL
    mf.conv_tol = SCF_CONV_TOL
    mf.max_cycle = 200

    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"SCF not converged for {FUNCTIONAL}/{basis}")

    hess = mf.Hessian().kernel()
    ha = harmonic_analysis(mol, hess)
    freqs = [float(x) for x in ha["freq_wavenumber"]]
    n_imag = int(sum(1 for f in freqs if f < 0))
    min_freq = float(min(freqs))

    zpe_ha = zpe_from_freqs_cm1(freqs)

    return {
        "e_tot_ha": float(e),
        "freqs_cm1": freqs,
        "n_imag": n_imag,
        "min_freq_cm1": min_freq,
        "zpe_ha": zpe_ha,
        "zpe_kcal_mol": zpe_ha * HARTREE_TO_KCAL_MOL,
        "pyscf_version": pyscf_version,
    }


def main():
    data = json.loads(SUMMARY_PATH.read_text())
    # map basis -> xyz
    basis_to_xyz = {
        "def2-svp": WORKDIR / f"{NAME}_prod_def2-svp_opt.xyz",
        "def2-tzvp": WORKDIR / f"{NAME}_prod_def2-tzvp_opt.xyz",
    }

    for r in data["results"]:
        basis = r["basis"]
        out = run_for_basis(basis, basis_to_xyz[basis])
        # update
        r["e_tot_ha"] = out["e_tot_ha"]
        r["freqs_cm1"] = out["freqs_cm1"]
        r["n_imag"] = out["n_imag"]
        r["min_freq_cm1"] = out["min_freq_cm1"]
        r["zpe_ha"] = out["zpe_ha"]
        r["zpe_kcal_mol"] = out["zpe_kcal_mol"]
        r["pyscf_version"] = out["pyscf_version"]

    SUMMARY_PATH.write_text(json.dumps(data, indent=2))
    print("Updated", SUMMARY_PATH)


if __name__ == "__main__":
    main()
