#!/usr/bin/env python
"""SA-CASSCF(6,5)/def2-TZVP + NEVPT2 for formaldehyde (H2CO) with NOON diagnostics.

Goal (per user request):
- Use an expanded active space CAS(6e,5o) intended to include sigma_CO, pi_CO, n_O, pi*_CO, sigma*_CO
- Run SA-CASSCF with 6 roots (S0 + 5)
- Run a separate multi-root CASCI at SA orbitals
- Run NEVPT2 per root on the CASCI reference
- Report per-root NOONs (active-space natural occupations) and PT2 correction magnitudes

Writes: results_cas65.json
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

from pyscf import gto, scf, mcscf, mrpt

HARTREE_TO_EV = 27.211386245988

WORKDIR = Path("research/formaldehyde_excited_state")
XYZ_PATH = WORKDIR / "formaldehyde_pbe0_def2-tzvp_opt.xyz"
BASIS = "def2-tzvp"
CHARGE = 0
SPIN = 0

NROOTS = 6
NELECAS = 6
NCAS = 5


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


def active_noons_for_root(ci_solver, ci_vec, ncas: int, nelecas: int) -> List[float]:
    """Return NOONs (natural occupation numbers) in active space for a given root."""
    # spin0 formaldehyde -> use direct_spin1
    rdm1 = ci_solver.make_rdm1(ci_vec, ncas, nelecas)
    # diagonalize 1-RDM to get NOONs
    w = np.linalg.eigvalsh(rdm1)
    w = np.array(w, dtype=float)
    w = np.sort(w)[::-1]
    return [float(x) for x in w]


def main() -> None:
    WORKDIR.mkdir(parents=True, exist_ok=True)
    if not XYZ_PATH.exists():
        raise FileNotFoundError(f"Missing optimized geometry: {XYZ_PATH}")

    atom = xyz_to_atom_str(XYZ_PATH.read_text())
    mol = gto.M(
        atom=atom,
        unit="Angstrom",
        basis=BASIS,
        charge=CHARGE,
        spin=SPIN,
        symmetry=False,
        verbose=4,
    )

    # RHF reference (gate)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("Gate fail: RHF did not converge")

    # --- Choose active orbitals via AVAS on C/O valence AOs ---
    # This aims to capture sigma/pi bonding/antibonding + O lone pair.
    # We include C and O valence shells (2s,2p). Hydrogens contribute little.
    from pyscf.mcscf import avas

    ao_labels = [
        "C 2s",
        "C 2p",
        "O 2s",
        "O 2p",
    ]

    # AVAS returns (ncas, nelecas, mo_coeff). We then enforce requested CAS size.
    ncas_avas, nelecas_avas, mo_avas = avas.avas(mf, ao_labels)
    # Gate: AVAS must at least provide enough orbitals to pick from
    if ncas_avas < NCAS:
        raise RuntimeError(f"AVAS suggested ncas={ncas_avas}, smaller than requested NCAS={NCAS}")
    if nelecas_avas < NELECAS:
        raise RuntimeError(f"AVAS suggested nelecas={nelecas_avas}, smaller than requested NELECAS={NELECAS}")

    # --- (A) SA-CASSCF orbital optimization ---
    mc = mcscf.CASSCF(mf, NCAS, NELECAS)
    mc.fcisolver.spin = 0
    mc.max_cycle_macro = 200
    mc.conv_tol = 1e-8

    mc.fcisolver.nroots = NROOTS
    weights = [1.0 / NROOTS] * NROOTS
    mc = mc.state_average_(weights)

    mc.kernel(mo_avas)

    sa_e_states = [float(x) for x in mc.e_states]

    # --- (B) CASCI multi-root at SA orbitals ---
    ci = mcscf.CASCI(mf, NCAS, NELECAS)
    ci.fcisolver.spin = 0
    ci.fcisolver.nroots = NROOTS
    e_tot, e_cas, ci_vecs, mo_ci, mo_energy = ci.kernel(mc.mo_coeff)
    e_tot = [float(x) for x in np.array(e_tot).ravel()]

    # Per-root NOONs in active space
    noons = []
    for r in range(NROOTS):
        noons.append(active_noons_for_root(ci.fcisolver, ci_vecs[r], NCAS, NELECAS))

    # --- (C) NEVPT2 per root on CASCI ---
    corr_ha = []
    e_nevpt2 = []
    errors = []
    for r in range(NROOTS):
        try:
            nev = mrpt.NEVPT(ci)
            nev.root = r
            corr = float(nev.kernel())
            corr_ha.append(corr)
            e_nevpt2.append(float(e_tot[r] + corr))
            errors.append("")
        except Exception as ex:
            corr_ha.append(float("nan"))
            e_nevpt2.append(float("nan"))
            errors.append(repr(ex))

    # excitation energies
    e0 = e_tot[0]
    casci_exc_ev = [float((x - e0) * HARTREE_TO_EV) for x in e_tot]

    e0n = e_nevpt2[0]
    nevpt2_exc_ev = [float((x - e0n) * HARTREE_TO_EV) for x in e_nevpt2]

    out: Dict[str, Any] = {
        "molecule": "formaldehyde",
        "basis": BASIS,
        "active_space": {
            "nelecas": NELECAS,
            "ncas": NCAS,
            "selection": {
                "method": "AVAS",
                "ao_labels": ao_labels,
                "avas_ncas": int(ncas_avas),
                "avas_nelecas": int(nelecas_avas),
            },
        },
        "sa_casscf": {
            "method": f"SA-CASSCF({NELECAS},{NCAS})/{BASIS}",
            "nroots": NROOTS,
            "weights": weights,
            "e_states_ha": sa_e_states,
            "excitation_ev": [float((x - sa_e_states[0]) * HARTREE_TO_EV) for x in sa_e_states],
        },
        "casci": {
            "method": f"CASCI({NELECAS},{NCAS})/{BASIS} @ SA-CASSCF orbitals",
            "nroots": NROOTS,
            "e_states_ha": e_tot,
            "excitation_ev": casci_exc_ev,
            "noons": noons,
        },
        "nevpt2": {
            "method": f"NEVPT2//CASCI({NELECAS},{NCAS})/{BASIS} @ SA-CASSCF orbitals",
            "corr_ha": corr_ha,
            "corr_ev": [float(x * HARTREE_TO_EV) for x in corr_ha],
            "e_states_ha": e_nevpt2,
            "excitation_ev": nevpt2_exc_ev,
            "errors": errors,
        },
    }

    (WORKDIR / "results_cas65.json").write_text(json.dumps(out, indent=2))

    # Print a compact table
    print(f"=== SA-CASSCF({NELECAS},{NCAS}) + NEVPT2 for H2CO | {BASIS} ===")
    print("Root  CASCI_eV   NEVPT2_eV   corr(Ha)   corr(eV)   NOONs")
    for r in range(NROOTS):
        print(
            f"{r:>4d}  {casci_exc_ev[r]:9.4f}  {nevpt2_exc_ev[r]:10.4f}  {corr_ha[r]:9.6f}  {corr_ha[r]*HARTREE_TO_EV:9.3f}  {noons[r]}"
        )


if __name__ == "__main__":
    main()
