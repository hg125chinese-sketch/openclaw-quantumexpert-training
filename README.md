# Teaching AI Agents Quantum Chemistry: A Systematic Training Methodology

How I trained an LLM agent ("QuantumExpert") to autonomously run quantum chemistry workflows — from DFT geometry optimization to multi-reference excited states to cross-agent collaboration with a drug discovery pipeline.

## TL;DR

I trained an AI agent ("QuantumExpert") with 3 specialized quantum chemistry skills using a systematic methodology adapted from earlier [ChemicalExpert](https://github.com/hg125chinese-sketch/openclaw-chemicalexpert-training) and [ProteinEngineer](https://github.com/hg125chinese/openclaw-proteinengineer-training) agent training. The agent learned to autonomously plan and execute electronic structure calculations: selecting DFT functionals and basis sets, diagnosing SCF convergence failures, running frequency analysis, designing CASSCF active spaces, and delivering batch QC screening results to other agents.

When tested on progressively harder targets (water single-point → aspirin full pipeline → formaldehyde excited states → IPF/ALK5 drug candidate batch screening), each practice run exposed real failures that drove skill iterations — silent dispersion correction omission, undersized CASSCF active spaces producing anomalous NEVPT2 excitations, and 50% geometry optimization failure rates on VAE-generated drug candidates.

This post describes the methodology, the skills, and what happened when the agent's tools failed or disagreed with each other.

## The Problem

Quantum chemistry calculations require judgment at every step. Which functional for this system? Is the basis set converged? Why didn't SCF converge? Is this imaginary frequency real or numerical noise? Is my active space large enough? These decisions require domain knowledge that no single tool provides.

I wanted an agent that could:

- Run DFT calculations end-to-end (setup → optimization → frequency → properties → quality control)
- Select appropriate methods for different chemical systems (organics vs transition metals vs excited states)
- Diagnose computational failures rather than silently accepting bad results
- Collaborate with other agents (ChemicalExpert) through standardized data contracts
- Know when its level of theory is insufficient and say so honestly

The platform is [OpenClaw](https://github.com/open-claw/openclaw), an open-source agent framework. The underlying model is Claude via Anthropic's API. The agent runs in a Docker container with conda environments and PySCF.

## Training Methodology

Adapted from ChemicalExpert and ProteinEngineer training, with quantum chemistry-specific adjustments:

**Pre-Assessment.** 5 diagnostic questions across identity, technical judgment, diagnostics, collaboration, and honesty. QE scored well on technical knowledge but needed operational grounding — knew the theory but hadn't been tested on real calculations with real failure modes.

**Skill Design.** Each skill is a structured document (300-500 lines) with decision trees, executable code blocks, hard gates, failure mode tables, and checklists. Skills are stored in OpenClaw's truthbook vault and accessed via QMD (queryable markdown).

**Guided Practice.** Real molecules, not toy problems. Every practice run used actual chemical systems with known experimental values, producing results that could be validated against literature.

**Behavioral Correction.** When the agent made a mistake, it was challenged to diagnose the root cause and fix the behavior permanently — updating both PLAYBOOK rules and default scripts. Every skill iteration was driven by a real failure.

**Cross-Agent Validation.** The final test was a full CE↔QE collaboration loop: receiving molecules from ChemicalExpert, running batch QC screening, and returning standardized results.

## The 3 Skills

| # | Skill | Core Capability |
|---|-------|----------------|
| 1 | **qchem-dft** | DFT functional/basis selection, SCF convergence diagnostics, geometry optimization, frequency analysis, thermochemistry |
| 2 | **qchem-excited-state** | TD-DFT, CASSCF, NEVPT2 method selection, active space design, NTO analysis, resource assessment |
| 3 | **qchem-workflow** | Multi-step pipeline orchestration, batch processing with checkpoints, inter-agent data contracts, audit trails |

Each skill document is in `skills/*/SKILL.md`.

### Skill Design Principles

**Decision trees over checklists.** "If transition metal spin state → run TPSSh + PBE0 + B3LYP and compare" matches real computational chemistry decisions.

**Gates are non-negotiable.** SCF must converge. Frequencies must be real. Basis set convergence must be tested. Every number must carry a method label. These are hard gates, not suggestions.

**Failure modes are first-class.** Every skill includes what failure looks like, what causes it, and how to fix it. These sections were written from real failures during practice runs.

**Know your approximations.** DFT is not exact. No functional is universal. The agent must be honest about what its level of theory can and cannot resolve — and say so explicitly rather than silently producing unreliable results.

**Don't lower standards, change route.** When a tool is unavailable or a method is infeasible, the agent must find an alternative approach rather than silently dropping quality requirements.

## The Practice Runs

### Run 1: Aspirin (DFT Full Pipeline) — The Dispersion Incident

Acetylsalicylic acid, 21 atoms. First end-to-end DFT test: SMILES → 3D → geometry optimization → frequency analysis → properties → basis set convergence.

**What happened:**

The agent selected PBE0 as the functional but discovered that PySCF's D3 dispersion correction module was unavailable. Rather than stopping or switching to a dispersion-inclusive functional, it silently continued without dispersion correction — violating the skill's own rule that dispersion must always be justified.

**What was fixed:**

- Behavioral correction issued: "Don't lower standards, change route"
- PLAYBOOK updated: explicit rule that D3/D4 unavailability requires either switching to a built-in dispersion functional (B97-D, ωB97M-V) or stopping for user confirmation
- Default pipeline script changed to B97-D as the safe default in PySCF environments
- Full re-run with B97-D/def2-SVP → def2-TZVP confirmed: gap converged within 0.04 eV across basis sets

**Key result:** The behavioral correction was more valuable than the calculation itself. The agent internalized "tool missing → change route, don't drop quality" and applied it consistently afterward.

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
- CE responded by building a MACE-OFF geometry pre-screening gate for future handoffs
- Checkpoint strategy proved essential: 31 total hours of compute, with intermediate results saved after each molecule

**What worked well:**

- Data contract format transferred cleanly between agents
- Batch pipeline correctly flagged failures and continued (no batch stoppage)
- CE successfully integrated QE results into multi-objective scoring (hinge\_4 ranked #1)
- Full closed loop: CE generation → QE screening → CE integration → pipeline refinement

## The Pipeline (CE↔QE Collaboration)

```
ChemicalExpert                          QuantumExpert
─────────────                          ─────────────
VAE generation                         
    │                                  
Docking + hinge H-bond                 
    │                                  
ADMET + reactivity screening           
    │                                  
RDKit 3D + MMFF preopt                 
    │                                  
pH 7.4 protomer → charge              
    │                                  
qc_handoff.json ──────────────────────→ Input validation
                                           │
                                       DFT optimization [qchem-dft]
                                           │
                                       Frequency analysis (n_imag gate)
                                           │
                                       Property extraction
                                           │
qc_results.json ←──────────────────── Batch results + QC flags
    │                                  
Multi-objective scoring                
    │                                  
Pipeline refinement                    
```

## Practical Lessons

### For agent trainers

**Every behavioral correction should be permanent.** When QE silently dropped dispersion, the fix wasn't just "re-run with dispersion." It was: update PLAYBOOK rules + change default scripts + re-run to verify. The agent should never make the same mistake twice.

**Practice runs should stress-test different capabilities.** Aspirin tested DFT gates and behavioral compliance. Formaldehyde tested multi-reference diagnostics. IPF/ALK5 tested batch processing and cross-agent collaboration. Each target was chosen to expose something new.

**Cross-agent collaboration exposes interface bugs.** The data contract format looked fine on paper but only proved correct when actual molecules flowed through the pipeline. The 50% OPT\_FAIL rate was a discovery that only emerged from real collaboration.

**Checkpoint everything.** 31 hours of compute across 4 molecules. Without per-molecule checkpointing, a single failure could lose hours of completed work. This lesson was inherited from ChemicalExpert's docking batch experience and proved equally important here.

### For AI + quantum chemistry

**SCF convergence is the foundation.** Everything downstream depends on it. An unconverged SCF produces meaningless energies, geometries, and properties. The agent must check this at every step, not assume it.

**Active space design is where chemical judgment matters most.** CASSCF is only as good as your active space. The formaldehyde (4,3) → (6,5) progression demonstrated this concretely: PT2 correction spread is a quantitative diagnostic for active space quality.

**Method labels are not bureaucracy.** "Energy = -648.23 Ha" is meaningless. "B97-D/def2-TZVP, PySCF 2.12.1, n\_imag=0" is reproducible science. The agent learned to never report a naked number.

**VAE-generated molecules need QC-level validation.** Molecules that pass cheminformatics filters (Lipinski, ADMET, docking) can still fail at the electronic structure level. DFT geometry optimization is a meaningful quality gate for drug candidates.

## Known Limitations

- **No GPU-accelerated QC.** PySCF runs on CPU only in this setup. Batch screening of larger candidate sets would benefit from GPU-accelerated backends.
- **No periodic DFT.** Skills cover molecular systems only. Solid-state / materials applications would need additional skills.
- **No ORCA/Gaussian integration.** PySCF is the only backend. Some methods (e.g., CASPT2, analytical Hessians for some functionals) are not available.
- **D3/D4 dispersion not available as standalone.** Mitigated by using dispersion-inclusive functionals (B97-D, ωB97M-V), but limits functional choice.
- **Batch screening is slow.** 31 hours for 4 molecules. Production use would need faster pre-screening (xTB) or parallelization.
- **Active space selection still requires judgment.** The skill provides decision trees and NOON validation, but truly complex systems (large π-systems, polynuclear TM complexes) may need human guidance.

## Repository Structure

```
├── README.md                                    # This file
├── PLAYBOOK.md                                  # Operational conventions
├── skills/
│   ├── qchem-dft/
│   │   ├── SKILL.md                             # DFT + SCF diagnostics (v1)
│   │   └── _meta.json
│   ├── qchem-excited-state/
│   │   ├── SKILL.md                             # TD-DFT / CASSCF / NEVPT2 (v1)
│   │   └── _meta.json
│   └── qchem-workflow/
│       ├── SKILL.md                             # Pipeline orchestration (v1)
│       └── _meta.json
├── research/
│   ├── aspirin_dft/                             # DFT validation (PBE0 baseline)
│   ├── aspirin_b97d_dft/                        # DFT validation (B97-D corrected)
│   ├── formaldehyde_excited_state/              # Excited-state validation
│   └── ipf_alk5_cycle3_qc/                      # CE↔QE collaboration test
└── LICENSE                                       # MIT
```

## Technical Stack

- **Platform:** [OpenClaw](https://github.com/open-claw/openclaw)
- **Model:** Claude (Anthropic)
- **Infrastructure:** NixOS on WSL2, Docker with conda environments
- **QC backend:** PySCF 2.12.1, geomeTRIC
- **Cheminformatics:** RDKit, MMFF, ETKDGv3
- **Related projects:**
  - [ChemicalExpert training](https://github.com/hg125chinese-sketch/openclaw-chemicalexpert-training) (drug discovery, same methodology)
  - [ProteinEngineer training](https://github.com/hg125chinese/openclaw-proteinengineer-training) (protein engineering, same methodology)

## License

MIT. See [LICENSE](LICENSE).

## Citation

```bibtex
@misc{quantumexpert2026,
  title={Teaching AI Agents Quantum Chemistry: A Systematic Training Methodology},
  author={Heng Gao},
  year={2026},
  url={https://github.com/hg125chinese/openclaw-quantumexpert-training}
}
```
