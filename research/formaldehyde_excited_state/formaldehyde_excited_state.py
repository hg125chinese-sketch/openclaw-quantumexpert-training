#!/usr/bin/env python
"""Excited-state analysis for formaldehyde (H2CO) following qchem-excited-state skill.

Requested deliverables:
1) TD-PBE0/def2-TZVP: lowest 5 singlet excited states (E, lambda, f)
2) NTO analysis for each state; classify n->pi* vs pi->pi*
3) SA-CASSCF(4,3)/def2-TZVP comparison (active space: n + pi + pi*)
4) NEVPT2 corrections for each state
5) Compare TD-DFT vs NEVPT2 excitation energies and analyze differences

All outputs are written under WORKDIR.

Notes:
- We compute vertical excitation energies at the optimized *ground-state* geometry.
- Formaldehyde is small; def2-TZVP is feasible.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Any, List, Tuple

import numpy as np

from pyscf import __version__ as pyscf_version
from pyscf import gto, dft, tddft, scf, mcscf, mrpt
from pyscf.geomopt.geometric_solver import optimize

HARTREE_TO_EV = 27.211386245988
EV_TO_NM = 1239.8419843320026

WORKDIR = Path("research/formaldehyde_excited_state")
NAME = "formaldehyde"

BASIS = "def2-tzvp"
XC = "pbe0"
GRID_LEVEL = 4
SCF_CONV_TOL = 1e-10
TD_CONV_TOL = 1e-6
TD_MAX_CYCLE = 2000
TD_MAX_SPACE = 200
N_STATES = 5  # report first 5 singlet excited states
N_STATES_COMPUTE = 10  # compute more roots to improve Davidson convergence; report first N_STATES

CHARGE = 0
SPIN = 0

# A reasonable starting geometry (Å)
# Roughly: C=O ~1.21 Å, C-H ~1.10 Å, H-C-H angle ~116 deg
ATOM_INIT = """
C  0.000000  0.000000  0.000000
O  0.000000  0.000000  1.210000
H  0.000000  0.943000 -0.543000
H  0.000000 -0.943000 -0.543000
""".strip()


@dataclass
class TDState:
    idx: int
    e_ev: float
    lam_nm: float
    f_osc: float
    nto_w0: float
    classification: str
    note: str


def build_mol(atom: str) -> gto.Mole:
    return gto.M(
        atom=atom,
        unit="Angstrom",
        basis=BASIS,
        charge=CHARGE,
        spin=SPIN,
        symmetry=False,
        verbose=4,
    )


def run_ground_opt(atom: str) -> gto.Mole:
    mol0 = build_mol(atom)
    mf0 = dft.RKS(mol0)
    mf0.xc = XC
    mf0.grids.level = GRID_LEVEL
    mf0.conv_tol = SCF_CONV_TOL
    mf0.max_cycle = 200

    mol_opt = optimize(mf0, maxsteps=200)

    # Gate: SCF converged at optimized geometry
    mf = dft.RKS(mol_opt)
    mf.xc = XC
    mf.grids.level = GRID_LEVEL
    mf.conv_tol = SCF_CONV_TOL
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("Gate fail: ground-state SCF not converged at optimized geometry")

    return mol_opt


def write_xyz(mol: gto.Mole, path: Path, comment: str) -> None:
    bohr_to_ang = 0.529177210903
    coords = mol.atom_coords() * bohr_to_ang
    lines = [str(mol.natm), comment]
    for i in range(mol.natm):
        sym = mol.atom_symbol(i)
        x, y, z = coords[i]
        lines.append(f"{sym:2s} {x: .8f} {y: .8f} {z: .8f}")
    path.write_text("\n".join(lines) + "\n")


def ao_atom_slices(mol: gto.Mole) -> List[Tuple[str, slice]]:
    """Return list of (atom_symbol, slice of AO indices)"""
    aoslices = mol.aoslice_by_atom()
    out = []
    for ia in range(mol.natm):
        p0, p1 = aoslices[ia, 2], aoslices[ia, 3]
        out.append((mol.atom_symbol(ia), slice(p0, p1)))
    return out


def orbital_atom_pop(mol: gto.Mole, orb_ao: np.ndarray) -> Dict[str, float]:
    """Qualitative per-atom population for an AO-basis orbital.

    For robust classification (n vs π) we want *non-negative* atom-resolved weights.
    A strict Mulliken gross population can yield small negative atomic contributions.

    Here we use a simple, always-nonnegative proxy:
      pop_A ∝ sum_{mu in A} |C_mu|^2

    This is not a rigorous population analysis; it is only used to distinguish
    O-localized lone-pair character from C/O delocalized π character.
    """
    v = np.asarray(orb_ao, dtype=float)
    pops = {}
    for sym, sl in ao_atom_slices(mol):
        pops[sym] = float(np.sum(v[sl] ** 2))
    tot = sum(pops.values())
    if tot > 1e-16:
        pops = {k: val / tot for k, val in pops.items()}
    return pops


def classify_state_from_nto(mol: gto.Mole, mf: dft.rks.RKS, td: tddft.TDDFT, state: int) -> Tuple[float, str, str]:
    """Return (w0, classification, note) for excited state index starting at 1.

    We compute NTOs via an explicit SVD of the transition density in the MO basis.

    Rationale: PySCF's `td.get_nto()` return format can vary (sometimes returns a
    full MO-rotation matrix), which is inconvenient for automated classification.

    Procedure (dominant pair):
    - Get TD amplitudes (X,Y) for the state
    - Form transition density T = X + Y (MO occ->vir block)
    - SVD: T = U S V^T
      * hole NTO (occ space)  = U[:,0]
      * particle NTO (vir space) = V[:,0]
    - Convert to AO: C_occ @ U[:,0], C_vir @ V[:,0]
    - Classify by atomic gross populations: O-localized hole => n→π*, else π→π*

    Returns:
      w0: dominant NTO weight (normalized S^2 fraction)
    """
    X, Y = td.xy[state - 1]
    X = np.asarray(X)
    Y = np.asarray(Y)
    T = X + Y

    # SVD in MO space
    U, s, Vt = np.linalg.svd(T, full_matrices=False)
    s2 = s**2
    w0 = float(s2[0] / s2.sum()) if s2.sum() > 0 else 0.0

    nocc = int(np.count_nonzero(mf.mo_occ > 0))
    C_occ = np.asarray(mf.mo_coeff[:, :nocc], dtype=float)
    C_vir = np.asarray(mf.mo_coeff[:, nocc:], dtype=float)

    hole_ao = C_occ @ U[:, 0]
    part_ao = C_vir @ Vt.T[:, 0]

    pop_h = orbital_atom_pop(mol, hole_ao)
    pop_p = orbital_atom_pop(mol, part_ao)

    o_h = pop_h.get("O", 0.0)
    c_h = pop_h.get("C", 0.0)

    if o_h > 0.70 and c_h < 0.25:
        cls = "n→π*"
        note = f"Dominant hole NTO localized on O (O={o_h:.2f}, C={c_h:.2f}); particle pops {pop_p}"
    else:
        cls = "π→π*"
        note = f"Dominant hole NTO not strongly O-localized (O={o_h:.2f}, C={c_h:.2f}); particle pops {pop_p}"

    return w0, cls, note


def run_tddft(mol: gto.Mole) -> Tuple[dft.rks.RKS, tddft.TDDFT, List[TDState]]:
    mf = dft.RKS(mol)
    mf.xc = XC
    mf.grids.level = GRID_LEVEL
    mf.conv_tol = SCF_CONV_TOL
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("Gate fail: ground-state SCF not converged before TD-DFT")

    td = tddft.TDDFT(mf)
    td.nstates = N_STATES_COMPUTE
    td.conv_tol = TD_CONV_TOL
    td.max_cycle = TD_MAX_CYCLE
    if hasattr(td, "max_space"):
        td.max_space = TD_MAX_SPACE
    # Provide explicit initial guess (can help Davidson converge)
    x0 = td.get_init_guess(mf, nstates=td.nstates)
    td.kernel(x0=x0)

    # Gate: all *reported* TD states must converge
    if hasattr(td, "converged"):
        conv = np.array(td.converged, dtype=bool)
        if conv.size >= N_STATES and (not np.all(conv[:N_STATES])):
            bad = (np.where(~conv[:N_STATES])[0] + 1).tolist()  # 1-based state indices
            raise RuntimeError(f"Gate fail: TD-DFT did not converge for reported states {bad}")

    # Oscillator strengths
    f = np.array(td.oscillator_strength(), dtype=float)

    states: List[TDState] = []
    for i in range(N_STATES):
        e_ev = float(td.e[i] * HARTREE_TO_EV)
        lam_nm = float(EV_TO_NM / e_ev) if e_ev > 0 else float("inf")
        f_osc = float(f[i])

        w0, cls, note = classify_state_from_nto(mol, mf, td, state=i + 1)

        states.append(
            TDState(
                idx=i + 1,
                e_ev=e_ev,
                lam_nm=lam_nm,
                f_osc=f_osc,
                nto_w0=w0,
                classification=cls,
                note=note,
            )
        )

    return mf, td, states


def run_sa_casscf_and_nevpt2(mol: gto.Mole, nroots: int = 6) -> Dict[str, Any]:
    """SA-CASSCF(4,3) orbital optimization + NEVPT2 corrections per root.

    IMPORTANT gate / PySCF constraint:
    - PySCF NEVPT2 cannot be run directly on a *state-averaged FCI solver*.
      We therefore:
      (A) run SA-CASSCF to obtain a balanced orbital set
      (B) run a separate multi-root CASCI (NOT state-averaged solver) on those orbitals
      (C) run NEVPT2 on that CASCI object for each root.

    Returns dict with:
      - SA-CASSCF root energies (from the SA-CASSCF itself; mostly for traceability)
      - CASCI root energies at the SA orbitals (used for NEVPT2 reference)
      - NEVPT2 corrected energies and excitation energies
    """
    # RHF reference (gate)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("Gate fail: RHF not converged for CASSCF reference")

    ncas = 3
    nelecas = 4

    # --- (A) SA-CASSCF orbital optimization ---
    mc = mcscf.CASSCF(mf, ncas, nelecas)
    mc.max_cycle_macro = 200
    mc.conv_tol = 1e-8
    mc.fcisolver.spin = 0  # singlets

    # Choose active orbitals explicitly for H2CO: pi (HOMO-1), n (HOMO), pi* (LUMO)
    nocc = mol.nelectron // 2
    cas_list = [nocc - 2, nocc - 1, nocc]  # 0-based indices in mo_coeff
    mo0 = mcscf.sort_mo(mc, mf.mo_coeff, cas_list)

    mc.fcisolver.nroots = nroots
    weights = [1.0 / nroots] * nroots
    mc = mc.state_average_(weights)
    mc.kernel(mo0)

    sa_e_roots = [float(x) for x in mc.e_states]

    # Orbitals to use for (B)
    mo_sa = mc.mo_coeff

    # --- (B) Multi-root CASCI at SA orbitals (NOT state-averaged solver) ---
    ci = mcscf.CASCI(mf, ncas, nelecas)
    ci.fcisolver.spin = 0
    ci.fcisolver.nroots = nroots
    e_tot, e_cas, ci_vecs, mo_ci, mo_energy = ci.kernel(mo_sa)

    e_roots = [float(x) for x in np.array(e_tot).ravel()]
    e0 = e_roots[0]
    casci_exc_ev = [float((er - e0) * HARTREE_TO_EV) for er in e_roots]

    # --- (C) NEVPT2 per root on CASCI object ---
    nevpt2_corr = []
    nevpt2_e = []
    nevpt2_err = []
    for r in range(nroots):
        try:
            nev = mrpt.NEVPT(ci)
            nev.root = r
            corr = float(nev.kernel())
            nevpt2_corr.append(corr)
            nevpt2_e.append(float(e_roots[r] + corr))
            nevpt2_err.append("")
        except Exception as ex:
            nevpt2_corr.append(float("nan"))
            nevpt2_e.append(float("nan"))
            nevpt2_err.append(repr(ex))

    nevpt2_exc_ev = [float((er - nevpt2_e[0]) * HARTREE_TO_EV) for er in nevpt2_e]

    return {
        "basis": BASIS,
        "casscf": {
            "method": f"SA-CASSCF({nelecas},{ncas})/{BASIS}",
            "nroots": nroots,
            "weights": weights,
            "e_states_ha": sa_e_roots,
            "excitation_ev": [float((er - sa_e_roots[0]) * HARTREE_TO_EV) for er in sa_e_roots],
        },
        "casci": {
            "method": f"CASCI({nelecas},{ncas})/{BASIS} @ SA-CASSCF orbitals",
            "nroots": nroots,
            "e_states_ha": e_roots,
            "excitation_ev": casci_exc_ev,
        },
        "nevpt2": {
            "method": f"NEVPT2//CASCI({nelecas},{ncas})/{BASIS} @ SA-CASSCF orbitals",
            "corr_ha": nevpt2_corr,
            "e_states_ha": nevpt2_e,
            "excitation_ev": nevpt2_exc_ev,
            "errors": nevpt2_err,
        },
        "active_space_note": "Active space fixed as (4e,3o) intended to span n + π + π*; selected as [HOMO-1, HOMO, LUMO] from RHF canonical orbitals.",
    }


def main() -> None:
    WORKDIR.mkdir(parents=True, exist_ok=True)

    # Phase A: ground-state optimization (gate upstream)
    mol_opt = run_ground_opt(ATOM_INIT)
    xyz_path = WORKDIR / f"{NAME}_{XC}_{BASIS}_opt.xyz"
    write_xyz(
        mol_opt,
        xyz_path,
        f"Optimized {NAME} at {XC}/{BASIS}, PySCF {pyscf_version}",
    )

    # Phase B: TD-DFT at optimized geometry
    mf, td, td_states = run_tddft(mol_opt)

    # Phase C: SA-CASSCF + NEVPT2
    mr = run_sa_casscf_and_nevpt2(mol_opt, nroots=N_STATES + 1)  # ground + 5

    # Assemble
    out = {
        "molecule": NAME,
        "atom_opt_xyz": str(xyz_path),
        "charge": CHARGE,
        "spin": SPIN,
        "software": {"pyscf": pyscf_version},
        "ground_state": {
            "xc": XC,
            "basis": BASIS,
            "grid_level": GRID_LEVEL,
            "scf_conv_tol": SCF_CONV_TOL,
        },
        "tddft": {
            "method": f"TD-{XC}/{BASIS}",
            "nstates": N_STATES,
            "td_conv_tol": TD_CONV_TOL,
            "states": [asdict(s) for s in td_states],
        },
        "multiref": mr,
    }

    (WORKDIR / "results.json").write_text(json.dumps(out, indent=2))

    # Human-readable summary
    lines = []
    lines.append(f"=== {NAME}: excited-state summary (vertical) ===")
    lines.append(f"Ground geom: {XC}/{BASIS}, PySCF {pyscf_version}")
    lines.append("")
    lines.append(f"TD-DFT: TD-{XC}/{BASIS} (nstates={N_STATES})")
    for s in td_states:
        lines.append(
            f"  S{s.idx}: E={s.e_ev:.4f} eV  (λ={s.lam_nm:.1f} nm)  f={s.f_osc:.4f}  NTO w0={s.nto_w0:.3f}  {s.classification}"
        )
    lines.append("")
    lines.append(mr["casscf"]["method"])
    for i, eev in enumerate(mr["casscf"]["excitation_ev"][1:], start=1):
        lines.append(f"  Root {i}: ΔE={eev:.4f} eV")
    lines.append("")
    lines.append(mr["nevpt2"]["method"])
    for i, eev in enumerate(mr["nevpt2"]["excitation_ev"][1:], start=1):
        lines.append(f"  Root {i}: ΔE={eev:.4f} eV")

    (WORKDIR / "summary.txt").write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    print(f"\nWrote: {WORKDIR}/results.json and summary.txt")


if __name__ == "__main__":
    main()
