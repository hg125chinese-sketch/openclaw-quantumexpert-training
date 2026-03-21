# Teaching AI Agents Quantum Chemistry: A Systematic Training Methodology

How I trained an LLM agent ("QuantumExpert") to autonomously run quantum chemistry workflows — from DFT geometry optimization to multi-reference excited states to cross-agent collaboration with a drug discovery pipeline.

## TL;DR

I trained an AI agent ("QuantumExpert") with 4 specialized quantum chemistry skills using a systematic methodology adapted from earlier [ChemicalExpert](https://github.com/hg125chinese-sketch/openclaw-chemicalexpert-training) and [ProteinEngineer](https://github.com/hg125chinese-sketch/openclaw-proteinengineer-training) agent training. The agent learned to autonomously plan and execute electronic structure calculations: selecting DFT functionals and basis sets, diagnosing SCF convergence failures, running frequency analysis, designing CASSCF active spaces, and delivering batch QC screening results to other agents.

Over seven DMTA cycles, an analog exploration campaign, N-N-free de-risking, and successor optimization on an IPF/ALK5 drug discovery target, the CE↔QE collaboration processed 20 molecules through DFT screening — achieving 100% pass rate across the entire DiffSBDD era (14/14). The agent evolved from a B97-D workaround (D3 unavailable) to B3LYP-D3(BJ) as the verified default, reducing DFT wall time from 7-10 hours to under 2 hours per molecule. The final molecule (NNF05_S05) achieved the project-wide highest HOMO-LUMO gap (4.97 eV) and lowest dipole moment (1.54 D) — the cleanest electronic structure profile in the entire campaign.

All training artifacts, skills, and case studies are open-source.

## The Problem

Quantum chemistry calculations require judgment at every step. Which functional for this system? Is the basis set converged? Why didn't SCF converge? Is this imaginary frequency real or numerical noise? Is my active space large enough? These decisions require domain knowledge that no single tool provides.

I wanted an agent that could:

- Run DFT calculations end-to-end (setup → optimization → frequency → properties → quality control)
- Select appropriate methods for different chemical systems (organics vs transition metals vs excited states)
- Diagnose computational failures rather than silently accepting bad results
- Collaborate with other agents (ChemicalExpert) through standardized data contracts
- Know when its level of theory is insufficient and say so honestly

The platform is [OpenClaw](https://github.com/open-claw/openclaw), an open-source agent framework. The underlying model is GPT-5.4 via OpenAI's Codex API. The agent runs in a Docker container with conda environments, PySCF, and xTB.

## Training Methodology

Adapted from ChemicalExpert and ProteinEngineer training, with quantum chemistry-specific adjustments:

**Pre-Assessment.** 5 diagnostic questions across identity, technical judgment, diagnostics, collaboration, and honesty. QE scored well on technical knowledge but needed operational grounding — knew the theory but hadn't been tested on real calculations with real failure modes.

**Skill Design.** Each skill is a structured document (200-500 lines) with decision trees, executable code blocks, hard gates, failure mode tables, and checklists. Skills are stored in OpenClaw's truthbook vault and accessed via QMD (queryable markdown).

**Guided Practice.** Real molecules, not toy problems. Every practice run used actual chemical systems with known experimental values, producing results that could be validated against literature.

**Behavioral Correction.** When the agent made a mistake, it was challenged to diagnose the root cause and fix the behavior permanently — updating both PLAYBOOK rules and default scripts. Every skill iteration was driven by a real failure.

**Cross-Agent Validation.** The final test was a full CE↔QE collaboration loop: receiving molecules from ChemicalExpert, running batch QC screening, and returning standardized results.

## The 4 Skills

| # | Skill | Core Capability |
|---|-------|----------------|
| 1 | **qchem-dft** | DFT functional/basis selection, SCF convergence diagnostics, geometry optimization, frequency analysis, thermochemistry. B3LYP-D3(BJ) as verified default via pyscf-dispersion. |
| 2 | **qchem-excited-state** | TD-DFT, ADC(2) via adcc (verified), multi-reference escalation (SA-CASSCF + NEVPT2/DF-NEVPT2), active space design with NOON validation, NTO analysis, resource assessment |
| 3 | **qchem-workflow** | Multi-step pipeline orchestration, batch processing with checkpoints, inter-agent data contracts (CE→QE and QE→CE JSON schemas), DiffSBDD 3D geometry handoff support, audit trails |
| 4 | **qchem-prescreen-xtb** | Fast GFN2-xTB prescreening (geometry sanity check before DFT), fail-fast heuristic with known limitation: xTB PASS does not guarantee DFT PASS |

Each skill document is in `skills/*/SKILL.md`.

### Skill Design Principles

**Decision trees over checklists.** "If transition metal spin state → run TPSSh + PBE0 + B3LYP and compare" matches real computational chemistry decisions.

**Gates are non-negotiable.** SCF must converge. Frequencies must be real. Basis set convergence must be tested. Every number must carry a method label. These are hard gates, not suggestions.

**Failure modes are first-class.** Every skill includes what failure looks like, what causes it, and how to fix it. These sections were written from real failures during practice runs.

**Know your approximations.** DFT is not exact. No functional is universal. The agent must be honest about what its level of theory can and cannot resolve — and say so explicitly rather than silently producing unreliable results.

**Don't lower standards, change route.** When a tool is unavailable or a method is infeasible, the agent must find an alternative approach rather than silently dropping quality requirements. This principle was established after the dispersion incident (Run 1) and permanently encoded in QE's PLAYBOOK.

## The Practice Runs

### Run 1: Aspirin (DFT Full Pipeline) — The Dispersion Incident

Acetylsalicylic acid, 21 atoms. First end-to-end DFT test: SMILES → 3D → geometry optimization → frequency analysis → properties → basis set convergence.

**What happened:**

The agent selected PBE0 as the functional but discovered that PySCF's D3 dispersion correction module was unavailable. Rather than stopping or switching to a dispersion-inclusive functional, it silently continued without dispersion correction — violating the skill's own rule that dispersion must always be justified.

**What was fixed:**

- Behavioral correction issued: "Don't lower standards, change route"
- PLAYBOOK updated: explicit rule that D3/D4 unavailability requires either switching to a built-in dispersion functional (B97-D, ωB97M-V) or stopping for user confirmation
- Default pipeline script changed to B97-D as the safe default
- Full re-run with B97-D/def2-SVP → def2-TZVP confirmed: gap converged within 0.04 eV across basis sets

**Post-training update:** pyscf-dispersion was later confirmed available in the environment, restoring B3LYP-D3(BJ) as the verified default. B97-D was demoted to fallback. DFT wall time per molecule dropped from 7-10 hours to under 2 hours.

### Run 2: Formaldehyde (Excited States) — Active Space Too Small

H₂CO, the textbook excited-state molecule. TD-DFT + CASSCF/NEVPT2 comparison with experimental S₁ (n→π*) ≈ 4.0 eV.

**What happened:**

- TD-PBE0/def2-TZVP gave S₁ = 4.043 eV — excellent agreement with experiment (+0.04 eV)
- SA-CASSCF(4,3)/def2-TZVP + NEVPT2 produced Root1 = 0.97 eV — absurdly low

**The diagnosis:**

The agent correctly identified the problem without prompting: the (4,3) active space (n + π + π\*) was too small. PT2 corrections were wildly state-dependent (spread of 2.14 eV across roots), causing unphysical root reordering. The skill's own active space guide warned that carbonyls "may need σ\_CO."

**What was fixed:**

- Expanded to CASSCF(6,5) including σ\_CO and σ\*\_CO
- PT2 correction spread dropped from 2.14 eV to 0.73 eV
- NEVPT2 Root2 = 4.10 eV, now consistent with experiment and TD-DFT
- NOON analysis confirmed all 5 orbitals were genuinely active

**Key result:** The diagnostic ability — recognizing that anomalous PT2 corrections signal an inadequate active space — was exactly the judgment the skill was designed to build.

### Run 3: IPF/ALK5 CE↔QE Collaboration — Batch Screening Reality Check

4 drug candidate molecules from ChemicalExpert's VAE generative model (IPF/ALK5 kinase inhibitor project, Cycle 3). First real cross-agent collaboration test.

**What happened:**

- CE preprocessed Top5 candidates: RDKit standardization, ETKDGv3 embedding, MMFF optimization, pH 7.4 charge determination
- CE packaged molecules in the qchem-workflow data contract format (JSON with mol\_id, SMILES, charge, spin)
- QE ran batch pipeline: B97-D/def2-SVP geometry optimization → frequency analysis → property extraction
- Results: 2/4 PASS, 2/4 OPT\_FAIL (geometry optimization did not converge)
- 50% failure rate exceeded the 20% warning threshold defined in qchem-workflow

**What was learned:**

- VAE-generated molecules can be chemically valid but geometrically pathological at the DFT level
- MMFF pre-optimization convergence does not guarantee DFT convergence
- xTB prescreening was validated but shown to not predict DFT failure (4/4 xTB PASS, 2/4 DFT OPT_FAIL)
- CE responded by building a MACE-OFF geometry pre-screening gate for future handoffs
- Checkpoint strategy proved essential: 31 total hours of compute, with intermediate results saved after each molecule

## The Pipeline (CE↔QE Collaboration)

```
ChemicalExpert                          QuantumExpert
─────────────                          ─────────────
DiffSBDD generation
    │
PoseBusters QC
    │
Safety screen (3-layer)
    │
Docking (Vina + multi-seed)
    │
GNINA CNN rescoring
    │
Interaction analysis (ProLIF)
    │
Panel selection (conflict-aware)
    │
RDKit 3D + MMFF preopt
    │
MACE geometry prescreen
    │
qc_handoff.json ──────────────────────→ Input validation
                                           │
                                       (optional) xTB prescreen
                                           │
                                       DFT optimization [qchem-dft]
                                       B3LYP-D3(BJ)/def2-SVP
                                           │
                                       Frequency analysis (n_imag gate)
                                           │
                                       Property extraction
                                           │
qc_results.json ←──────────────────── Batch results + QC flags
    │
Multi-objective scoring
    │
Retrosynthesis + analog exploration
    │
N-N-free de-risking + successor optimization
```

## CE↔QE Collaboration Statistics

Over seven DMTA cycles, analog exploration, N-N-free de-risking, and successor optimization:

| Campaign | Method | Sent | PASS | Fail Rate |
|----------|--------|------|------|-----------|
| Cycle 3 | B97-D/def2-SVP | 4 | 2 | 50% |
| Cycle 4 | B97-D/def2-SVP | 2 | 1 | 50% |
| Cycle 5 | B97-D/def2-SVP | 3 | 3 | 0% |
| Cycle 6 | B97-D/def2-SVP | 2 | 2 | 0% |
| Cycle 7 | B3LYP-D3(BJ)/def2-SVP | 1 | 1 | 0% |
| Analogs (mol\_0021) | B3LYP-D3(BJ)/def2-SVP | 3 | 3 | 0% |
| NNF (N-N-free) | B3LYP-D3(BJ)/def2-SVP | 3 | 3 | 0% |
| NNF05 successors | B3LYP-D3(BJ)/def2-SVP | 2 | 2 | 0% |
| **Total** | — | **20** | **17** | **15%** |

DiffSBDD era (Cycles 5-7 + analogs + NNF + NNF05 successors): **14/14 = 100% PASS**.

### QE Electronic Structure Highlights

| Molecule | Campaign | gap (eV) | dipole (D) | Vina | Note |
|----------|----------|----------|------------|------|------|
| NNF05\_S05 | NNF05 successor | **4.97** | **1.54** | -9.64 | Project-wide best gap + lowest dipole |
| NNF07 | NNF | 4.96 | 5.14 | — | Second highest gap |
| A3\_01 | Analog | 4.79 | 1.62 | — | Highest score\_final (10.635) |
| NNF05\_S10 | NNF05 successor | 4.76 | 2.67 | — | Strongest Boltz-2 among successors (0.591) |
| NNF02 | NNF | 4.61 | 3.95 | **-10.71** | Project-wide best Vina + N-N-free |

## Ecosystem Deep Research (2026)

After the initial training and collaboration cycles, QE conducted a systematic survey of the PySCF ecosystem (2024-2026) to identify upgrade opportunities:

| Discovery | Impact | Status |
|-----------|--------|--------|
| **pyscf-dispersion available** | B3LYP-D3(BJ) restored as default; B97-D demoted to fallback | ✅ Integrated |
| **adcc/ADC(2) works with PySCF 2.12.1** | Fills gap between TD-DFT and CASSCF as intermediate method | ✅ Verified, in skill |
| **DF-NEVPT2 in PySCF 2.12.0+** | ~1.4x speedup for multi-reference excited states | ✅ Verified, in skill |
| **GPU4PySCF** | 2.7x speedup on single-point (direct SCF); DF mode blocked by CuPy 14 | 📋 Pending CuPy fix |
| **ANI-2x / AIMNet2** | Multi-conformer prescreening; cannot predict DFT failure | 📋 Survey complete |
| **MC-PDFT (pyscf-forge)** | CASSCF+DFT hybrid; build failed locally (BLAS missing) | 📋 Pending |

## Practical Lessons

### For agent trainers

**Every behavioral correction should be permanent.** When QE silently dropped dispersion, the fix wasn't just "re-run with dispersion." It was: update PLAYBOOK rules + change default scripts + re-run to verify. The agent should never make the same mistake twice.

**Practice runs should stress-test different capabilities.** Aspirin tested DFT gates and behavioral compliance. Formaldehyde tested multi-reference diagnostics. IPF/ALK5 tested batch processing and cross-agent collaboration. Each target was chosen to expose something new.

**Cross-agent collaboration exposes interface bugs.** The data contract format looked fine on paper but only proved correct when actual molecules flowed through the pipeline. The 50% OPT\_FAIL rate was a discovery that only emerged from real collaboration.

**Checkpoint everything.** 31 hours of compute across 4 molecules. Without per-molecule checkpointing, a single failure could lose hours of completed work.

### For AI + quantum chemistry

**SCF convergence is the foundation.** Everything downstream depends on it. An unconverged SCF produces meaningless energies, geometries, and properties. The agent must check this at every step, not assume it.

**Active space design is where chemical judgment matters most.** CASSCF is only as good as your active space. The formaldehyde (4,3) → (6,5) progression demonstrated this concretely: PT2 correction spread is a quantitative diagnostic for active space quality.

**Method labels are not bureaucracy.** "Energy = -648.23 Ha" is meaningless. "B3LYP-D3(BJ)/def2-SVP, PySCF 2.12.1, n\_imag=0" is reproducible science.

**VAE-generated molecules need QC-level validation.** Molecules that pass cheminformatics filters (Lipinski, ADMET, docking) can still fail at the electronic structure level. DFT geometry optimization is a meaningful quality gate for drug candidates.

**xTB prescreening has limits.** xTB PASS does not predict DFT PASS (validated: 4/4 xTB PASS, 2/4 DFT OPT\_FAIL). Its value is fail-fast for catastrophic geometries and cheap geometry preconditioning, not convergence prediction.

**DiffSBDD 3D coordinates are for design, not QC.** Pocket-conditioned diffusion generates better SMILES but its 3D coordinates carry extreme force-field strain (450-630 kcal/mol). RDKit re-embed is required before DFT.

**DiffSBDD molecules are electronically cleaner.** 14/14 DFT PASS in the DiffSBDD era vs 3/6 in the VAE era. Pocket-conditioned generation produces molecules that are not only better docking candidates but also more tractable for quantum chemistry.

## Known Limitations

- **GPU4PySCF density fitting blocked.** Direct SCF works (2.7x speedup) but DF mode fails due to CuPy 14 incompatibility. Hessian acceleration unavailable in direct mode.
- **No ORCA/Gaussian integration.** PySCF is the only backend. Some methods (e.g., CASPT2, some analytical Hessians) are not available.
- **xTB does not predict DFT failure.** Validated empirically; kept as geometry preconditioner only.
- **Active space selection still requires judgment.** The skill provides decision trees and NOON validation, but truly complex systems may need human guidance.
- **Batch screening is CPU-bound.** ~2 hours per molecule with B3LYP-D3(BJ) (improved from 7-10h with B97-D). Production use would benefit from GPU acceleration or parallelization.

## Repository Structure

```
├── README.md                                    # This file
├── PLAYBOOK.md                                  # Operational conventions
├── skills/
│   ├── qchem-dft/
│   │   ├── SKILL.md                             # DFT + SCF diagnostics
│   │   └── _meta.json
│   ├── qchem-excited-state/
│   │   ├── SKILL.md                             # TD-DFT / ADC(2) / CASSCF / NEVPT2
│   │   └── _meta.json
│   ├── qchem-workflow/
│   │   ├── SKILL.md                             # Pipeline orchestration
│   │   └── _meta.json
│   └── qchem-prescreen-xtb/
│       ├── SKILL.md                             # xTB prescreening
│       └── _meta.json
├── research/
│   ├── aspirin_dft/                             # DFT validation (PBE0 baseline)
│   ├── aspirin_b97d_dft/                        # DFT validation (B97-D corrected)
│   ├── formaldehyde_excited_state/              # Excited-state validation
│   └── ipf_alk5_cycle3_qc/                      # CE↔QE collaboration test
└── LICENSE                                       # MIT
```

## Technical Stack

- **Platform:** [OpenClaw](https://github.com/open-claw/openclaw) v2026.3.11
- **Model:** GPT-5.4 via OpenAI Codex API
- **Infrastructure:** NixOS on WSL2, Docker with conda environments, RTX 4080 Laptop GPU
- **QC backend:** PySCF 2.12.1, pyscf-dispersion 1.5.0, geomeTRIC 1.1, xTB 6.7.1, adcc 0.16.1
- **Cheminformatics:** RDKit 2025.03.6, MMFF, ETKDGv3
- **Related projects:**
  - [ChemicalExpert training](https://github.com/hg125chinese-sketch/openclaw-chemicalexpert-training) (drug discovery, 28 skills, 7 cycles + analog + NNF campaigns)
  - [ProteinEngineer training](https://github.com/hg125chinese-sketch/openclaw-proteinengineer-training) (protein engineering, 7 skills, 3 practice runs)
  - [Multi-agent methodology overview](https://github.com/hg125chinese-sketch/openclaw-multiagent-drug-discovery)

## License

MIT. See [LICENSE](LICENSE).

## Citation

```bibtex
@misc{quantumexpert2026,
  title={Teaching AI Agents Quantum Chemistry: A Systematic Training Methodology},
  author={Heng Gao},
  year={2026},
  url={https://github.com/hg125chinese-sketch/openclaw-quantumexpert-training}
}
```
