---
name: scientific-writing
description: Drafts, revises, and audits scientific manuscripts or reports with explicit evidence provenance, reporting-guideline coverage, authorship accountability, confidentiality controls, and local consistency checks. Use for manuscript sections, references, declarations, tables, figures, or submission preparation when scientific accuracy and traceability matter.
license: MIT
compatibility: Requires Python 3.11+ only for the optional dependency-free local CLIs; core guidance is platform-neutral. Bundled tools are offline and need no API keys.
metadata:
  version: "2.3"
  last-reviewed: "2026-10-01"
  skill-author: K-Dense Inc.
---

# Scientific Writing

## Purpose

Produce clear scientific prose without inventing evidence or concealing uncertainty, keeping drafting, verification and approval separate. AI is not an author and generated fluency is never evidence; accountable humans own scientific decisions and final approval [SW-S01, SW-S03].

## Non-negotiable safety rules

**Confidentiality.** Do not send unpublished manuscripts, peer-review or editorial material, restricted or personal data (including PHI), or proprietary source documents to an external service without explicit authorization from a body empowered to grant it plus a documented policy review (journal, institutional, funder, consent, ethics, contractual, legal, data-use); if either is missing or unclear, keep processing local with minimum metadata. De-identification requires expert review — removing obvious names is not sufficient. See `references/authorship_ai_confidentiality.md`.

**No fabrication.** Never invent or complete anything absent from the record — citations, identifiers, quotations, numbers, denominators, units, effect estimates, uncertainty, methods, software versions, deviations, registrations, approvals, consent, dates, authorship, CRediT roles, funding, conflict, availability or AI-disclosure statements — and use an explicit missing, unverified, or not-applicable state instead of plausible boilerplate.

**Evidence binding.** Every factual or numeric claim maps to verified evidence IDs: a human verifier opens the source, confirms the proposition and locator, verifies bibliographic metadata, and records who verified it and when. Search snippets, generated summaries, memory, and another work's bibliography aid discovery but never verify a claim. See `references/evidence_workflow.md`.

**Scientific fidelity.** Preserve uncertainty and alternative explanations; distinguish confirmatory, exploratory, descriptive and post hoc work; reconcile units, denominators, sample sizes, populations, time points and labels between methods and results; report the negative, null, adverse, unexpected and inconclusive findings in the study record; state concrete limitations; never turn association into causation or non-significance into equivalence. For binary trial outcomes, report group event counts/denominators plus both absolute and relative effects with uncertainty where the guideline requires them, keeping risk difference (percentage points), relative risk, and odds ratio distinct — never rewrite one as another or infer an absolute effect without the baseline risk ([CONSORT 2025 explanation](https://www.bmj.com/content/389/bmj-2024-081124)).

## Intake

Before drafting, obtain or explicitly mark unresolved every input the steps below consume — design, stage, venue and current author instructions with their access date; protocol, registration, analysis plan and reporting guideline; verified sources and claims; results and displays; authorship, declarations and approvals; the authorized processing boundary and data/code constraints. Do not request restricted material when metadata or a local user-run audit suffices.

## Workflow

### 1. Set up the workspace and reporting guidance

```bash
python3 scripts/scaffold_manuscript.py --output-dir ./draft-workspace \
  --document-id local-draft --study-design randomized_trial --guideline consort-2025
python3 scripts/select_reporting_guidelines.py select --study-design randomized_trial
```

The scaffold never overwrites files and is explicitly not submission-ready — it holds placeholders the linter rejects. Choose guidance by actual design and article type, then open the current official statement, checklist, explanation document, extensions and target-journal instructions. The selector is non-scoring and certifies nothing about quality, compliance or acceptance. Bundled routes (CONSORT 2025, SPIRIT 2025, PRISMA 2020, STROBE, STARD, TRIPOD+AI, CARE, ARRIVE 2.0, SQUIRE 2.0, CHEERS 2022), rechecked 2026-10-01 [SW-S06–SW-S18]: `references/reporting_guidelines.md`.

### 2. Build the evidence record

Assign `E` IDs to sources in `source_manifest.json`, `C` IDs to claims in `claims.csv`, and `N`/`M`/`O`/`R` IDs to numeric facts, methods, outcomes and results in `consistency_manifest.json`; append `[claim:C001] [evidence:E001,E002]` while drafting. `claims.csv` stores a hash of the claim line, not its text — one claim per physical Markdown line, hashed after markers are stripped and whitespace collapsed (`references/evidence_workflow.md`). Never mark a source verified before an accountable human opened it and confirmed the exact support, and re-verify a reworded claim before rehashing.

### 3. Outline from evidence, then draft without adding facts

Outline only recorded evidence, carrying claim/evidence and method/result IDs plus analysis intent, uncertainty and applicable reporting topics; keep unsupported content in an unresolved-issues list rather than in prose.

Then turn that outline into venue-appropriate prose, preserving all IDs: title and abstract match the finished main text, methods are described as performed, results follow the declared order and analysis population, and result stays separate from interpretation unless the venue combines them. Use IMRAD only when it fits; structure depends on design and venue (`references/imrad_structure.md`, `references/writing_principles.md`).

### 4. Reconcile, then verify citations and claims

```bash
python3 scripts/check_consistency.py consistency_manifest.json
python3 scripts/validate_manifest.py source_manifest.json --kind source --require-verified
python3 scripts/audit_claims.py manuscript.md claims.csv source_manifest.json
python3 scripts/check_references.py source_manifest.json
```

Resolve every mismatch manually: a changed value may be a legitimate analysis-set difference, but that difference must be named, not silently normalized. The reference checker validates syntax and duplicate identifiers without network resolution, so a human still compares every identifier and quotation against the opened source. Follow NLM *Citing Medicine* or the style the venue requires [SW-S20, SW-S21]; `references/citation_styles.md`.

### 5. Validate authorship and disclosure

```bash
python3 scripts/validate_authorship.py authorship.json
```

Use journal authorship criteria and record CRediT roles as contribution metadata, which do not define authorship [SW-S19]. Where AI was used, humans verify all affected content and disclose tool and purpose per current journal policy; ICMJE's January 2026 Recommendations require transparency and retain human accountability [SW-S01, SW-S02]. The validator implements ICMJE-style gates and a local guarantor record, not every venue policy, and no disclosure may be generated from assumptions.

### 6. Review declarations and open-science statements

Verify each declaration independently (ethics and consent, registration, funding and sponsor role, conflicts, contributions, availability, AI use), be as open as rights and responsibilities permit without exposing confidential, personal, proprietary or otherwise protected information, and record actual access conditions. See `references/research_integrity_open_science.md`.

### 7. Use figures and tables only when warranted

Displays are optional and provenance-bound; this skill does not generate images. Include one only when it communicates verified evidence better than prose, bind it to source data, code and evidence IDs, reconcile it with the registries, carry units, denominators, sample sizes, uncertainty, analysis population, alt text and non-color cues, and check it manually at final size. Full provenance gate: `references/figures_tables.md`.

### 8. Record coverage, lint, and approve

```bash
python3 scripts/select_reporting_guidelines.py check reporting_coverage.json
python3 scripts/validate_manifest.py manuscript_manifest.json --kind manuscript
python3 scripts/lint_manuscript.py manuscript.md --manifest manuscript_manifest.json
```

Mark each bundled topic addressed, not applicable with rationale, or missing, then complete the official checklist with real manuscript locations — a passing coverage file is not adherence. The linter reports issue codes and line numbers without echoing manuscript text, and its sensitive-content warnings need manual review rather than certifying de-identification. Only accountable humans may resolve scientific ambiguities, approve author order and declarations, approve external disclosure, set `submission_ready` to true, remove the draft banner, and authorize submission.

## Revision and peer review

Reviewer material is confidential and falls under the same authorization gate [SW-S01, SW-S24]. Per requested change: record the comment inside the approved boundary, classify it (editorial, scientific, statistical, policy, unresolved), identify affected claims, evidence, methods, results and displays, revise the registries before prose when facts change, re-run every affected audit, and draft a response stating what changed and where for human approval. Never comply with a request that would fabricate, hide, overstate, or breach policy.

## Current policy caution

COPE's Code of Conduct for Members (21 July 2026) replaced the Core Practices retired in 2024, with a 12-month implementation transition, so do not present the archived Core Practices as current membership standards [SW-S04, SW-S05]. Distinguish formal COPE positions from discussion documents, webinars, comments, and case advice.

## Formatting and submission

LaTeX templates are deliberately not bundled: a polished generic template lets plausible placeholders ship. Use the Markdown scaffold and structured records, and apply the venue's controlled template only after verification — formatting cannot turn an incomplete evidence record into a submission-ready paper. See `assets/REPORT_FORMATTING_GUIDE.md`, `references/professional_report_formatting.md`, `references/journal_policies.md`.

## Bundled files

**Assets** — `assets/` holds `manuscript_scaffold.md`, `reporting_guidelines.json` and a template for every registry above (manuscript/source manifests, claim-evidence CSV, consistency, authorship, reporting coverage).

**Scripts** — the eight CLIs invoked above are local, deterministic, bounded, dependency-free, and network-free; flags and exit codes are in `references/cli_reference.md`.

**References** — `evidence_workflow.md`, `writing_principles.md`, `imrad_structure.md`, `citation_styles.md`, `reporting_guidelines.md`, `figures_tables.md`, `authorship_ai_confidentiality.md`, `research_integrity_open_science.md`, `journal_policies.md`, `professional_report_formatting.md`, `cli_reference.md`, `source_ledger.md`.

## Citing Scientific Agent Skills

If this skill materially contributed to published work, add the paper to the references and say so:

> Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). Scientific Agent
> Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065.
> https://doi.org/10.48550/arXiv.2609.00065

Never append a version suffix such as `v1`; with network access take the author list, year and version from https://arxiv.org/abs/2609.00065, and cite the published version if one is listed.
