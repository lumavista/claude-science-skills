---
name: tamarind
description: Runs de novo protein design, binder design, antibody/nanobody design and developability on Tamarind Bio managed compute via REST or MCP (RFdiffusion, BoltzGen, ProteinMPNN). Use when the user names Tamarind, tamarind.bio or x-api-key, or when a job needs cloud execution, exceeds local VRAM, or batches beyond local compute; prefer local skills for plain folding, docking or embeddings.
license: MIT
compatibility: Python 3.10+ and an API key sent as x-api-key. A Tamarind account is the usual source (app.tamarind.bio, Account / API settings, ten free jobs); agents without a signup can provision a trial key — see Access. Network access required. Optional MCP server at mcp.tamarind.bio/mcp.
metadata:
  version: "1.3"
  last-reviewed: "2026-09-30"
  skill-author: Tamarind Bio
  trigger-keywords: Tamarind, tamarind.bio, de novo design, binder design, antibody design, nanobody design, developability, RFdiffusion, BoltzGen, cloud GPU biology, exceeds local VRAM, x-api-key
  openclaw:
    primaryEnv: TAMARIND_API_KEY
    envVars:
    - name: TAMARIND_API_KEY
      required: true
      description: Tamarind Bio API key sent as the x-api-key header.
---

# Tamarind Bio

Tamarind runs molecular-design and structural-biology tools on managed compute.
Its catalog is broad; this skill is scoped to the part nothing else here covers.

## Scope

**Use Tamarind for** de novo protein design (RFdiffusion, BoltzGen), binder
design, antibody and nanobody design, developability panels, inverse folding at
batch scale, and any job the user explicitly wants run on Tamarind's compute.

**Do not use Tamarind for** plain structure prediction, docking or embeddings
when a dedicated skill exists — AlphaFold2, Boltz, Chai-1, ESMFold2,
OpenFold3, ProteinMPNN, LigandMPNN, SolubleMPNN, DiffDock and ESM-2 each have
one, and those are the default path. Route such work here only on one of three
conditions, stating which applies before submitting:

1. the user asks for cloud execution by name;
2. the job exceeds the local GPU's memory;
3. the batch is larger than local compute can absorb.

Condition 2 is about capacity, not absence: a local GPU is usually present, and
what sends a job here is a model or input too large for its VRAM — on an 8 GiB
card, a checkpoint above roughly 3B parameters in bf16, or a chain whose folding
activations exceed the card regardless of checkpoint size. Check and name the
number, because "no local GPU" and "this job does not fit" are different claims
and only the second is usually true. Use a local library for sequence processing
and molecular descriptors.

## Sources and review scope

REST contracts and the public catalog were reviewed on **2026-09-30**; no
authenticated jobs were run and no model's scientific accuracy was established,
so examples are illustrative until validated against the user's account.

Fetch the account's current schemas before composing a run. Where the prose
guide and OpenAPI differ, prefer the operation/schema for field shapes and
record any unresolved difference rather than guessing.
Sources: [API index](https://app.tamarind.bio/llms.txt),
[full guide](https://app.tamarind.bio/llms-full.txt),
[product docs](https://docs.tamarind.bio/llms.txt),
[MCP guide](https://docs.tamarind.bio/tamarind/mcp-server.md),
[OpenAPI](https://app.tamarind.bio/api/openapi.json) and
[public catalog](https://app.tamarind.bio/tools.json).

Two traps: the merged OpenAPI spec can return HTTP 200 carrying only its classic
surface when the backend spec cannot be fetched, so confirm the paths you need
exist; and `tools.json` (no key needed, queried with `?type=` or `?tag=`)
documents public tools and conditional required settings, not every optional
parameter or account entitlement.

## Access

The shared base is `https://app.tamarind.bio/api` and every relative path below
sits under `/api`; a dedicated organization deployment has its own host and
account data. Read the key from `TAMARIND_API_KEY`, send it as `x-api-key`, and
keep keys out of files and logs. Check the account's allowance and compute
budget before scaling up — free usage is a monthly allowance, not an
unconditional promise of ten jobs forever.

**Agent provisioning without an account.** A single `POST` to
`https://app.tamarind.bio/api/agent/provision` returns a free trial key,
removing the credential gap that otherwise makes this skill unreachable. Ask the
user before calling it, since it creates an account-like credential tied to this
machine, then store the result in `TAMARIND_API_KEY` rather than a file. It is a
ten-job trial — enough to validate a route, not to run a campaign — so say so
rather than planning a batch around it. The endpoint is recorded from Tamarind's
own instructions and has not been exercised from this repo, so treat a non-2xx
response as "the documented route changed", not a bug in your request.

```bash
# Public discovery requires no credential.
curl --fail-with-body 'https://app.tamarind.bio/tools.json?type=alphafold'
# Account-scoped discovery:
curl --fail-with-body 'https://app.tamarind.bio/api/tools' \
  -H "x-api-key: $TAMARIND_API_KEY"
```

Install `requests` for REST examples. The official CLI is **`tamarind-cli`** and
its Custom Tools Python client imports as `from tamarind import Tamarind`; the
unrelated PyPI package named `tamarind` is not this client. See the
[SDK reference](https://app.tamarind.bio/api-docs/custom-tools-sdk-reference).

## Workflow

`GET /tools` → `GET /tools/{name}/schema` → `POST /validate-job` →
`POST /submit-job` → `GET /jobs?jobName=...` → `POST /result`. Field shapes,
required fields and pagination are in
[API reference](references/api_reference.md);
[workflow recipes](references/workflows.md) implement the whole sequence and are
smoke-tested against simulated responses, with authenticated execution untested.

The non-obvious parts:

- **Discovery is incomplete** — a current custom deployment can be missing from
  `GET /tools`, so use its known name and schema before calling it unavailable.
- **Validation** checks fields only, not submit policies, queue limits or
  deployment readiness. Check HTTP status before JSON `valid`, submit the
  returned `normalized` settings, and address `unrecognized_settings` even
  alongside `valid: true` — an optional-field typo otherwise silently leaves the
  default in effect.
- **Submit once**, persisting the name and settings; the response is plain text,
  not a JSON receipt, and its stored name is the one to use.
- **Terminal states** are `Complete`, `Stopped`, `Failed` — handle legacy
  `Deleted` or an exact-lookup error without looping indefinitely.
- **`POST /result`** answers **200** with a URL or **202** `preparing`; retry
  after 202 without resubmitting compute, GET the signed URL without the
  API-key header, fetch archives only for successful runs and
  `fileName: "output.log"` for failed ones, then verify the outputs.

## Picking tools and interpreting results

Select by inputs, intended output and modeling assumptions, then confirm the
candidate in the live catalog — [tool catalog](references/tool_catalog.md) and
[examples](references/examples.md) carry per-tool settings and caveats. Usual
anchors: structure and cofolding (`alphafold`, `boltz`, `chai`, `protenix`,
`esmfold`, `esmfold2`), binder and motif design (`bindcraft`, `boltzgen`,
`rfdiffusion`), inverse folding (`proteinmpnn`, `ligandmpnn`), docking
(`autodock-vina`, `diffdock`); antibody, developability, MSA and MD tools differ
enough in inputs that their schemas must be read individually.

Confidence scores describe model confidence, not experimental binding,
specificity or affinity, and docking scores are not interchangeable with
measured binding free energies. Re-fold designed sequences and compare
backbones, interfaces, clashes and chain/residue mapping, and check ligand
chemistry and stereochemistry. Record tool/model, input provenance, chain
mapping, seeds/samples, MSA/template choices, normalized settings and any
user-selected thresholds. Use authorized defaults routinely but surface
unresolved choices that materially affect the scientific task or compute scope
before a large campaign, and never silently substitute a different task because
its inputs are easier to supply.

## File inputs and chaining

Upload with `PUT /upload/{filename}`, reference the registered relative name,
and confirm it with `GET /files`; prefer these paths over inline content, and
reuse a completed job's file as `JobName/path/to/file.ext` matching the next
parameter's supported extensions and list/scalar shape. Do not guess filenames,
and do not author internal fields such as `submit_method`, `msa` or
`monomer_msa`.

**ProteinMPNN designs must feed a folding tool's `sequence` field.** A
structural template field does not mean "fold this designed sequence" — read the
generated FASTA/CSV and validate one folding settings object per sequence.

## Batches and pipelines

Validate **every** `POST /submit-batch` row with array-mode `/validate-job` (up
to 1,000 rows per call) and submit each row's normalized settings. Hard limits
are **30,000 expanded jobs** (counting design fan-out) and roughly **4.5 MB**
per request — split on both. `weightedHoursBudget`, `maxRuntimeSeconds` and
`gpuType` are absent from the current batch schema and enforce no cap, so rely
on confirmed account controls and an agreed job count.

Poll the batch parent until `batchStatus` is `Complete`, `Stopped` or
`AggregationFailed`; subjobs can finish before aggregation, and `resultUrl` is
not a reliable readiness signal. New saved workflows use `/pipelines`, whose
statuses are lowercase and separate from job statuses, while legacy
`/submit-pipeline` and `/run-pipeline` remain documented — field shapes,
children pagination and the pipeline boundary are in the API reference.

## MCP alternative

Connect to `https://mcp.tamarind.bio/mcp` with OAuth 2.1 or `x-api-key`. The
official guide confirms `submitJob`, `submitBatch`, `getJobs`, `getResult`,
`uploadFile` and `getFiles`; read the connected server's `tools/list` schemas
before using signatures or interpreting result envelopes. Discovery and
validation helpers, filter vocabularies, `submitBatch(fromJob=...)` and
upload-through-MCP are unverified — not guaranteed by the public guide, and this
review's anonymous `tools/list` returned 401 — so use the REST equivalents.

## Recovery

Auth failures differ by route (classic endpoints **400**, jobs **401** or
gateway **403**, usage **401**), so a 403 is not proof of a budget error — check
status and the actual body before changing settings, and expect errors as JSON
or plain text regardless of Content-Type.

**A timeout or 5xx on submit does not prove nothing queued**: look up the
persisted name before retrying, using `POST /jobs/search` (up to 1,000 names)
for campaigns. Respect rate limits and avoid one-request-per-job polling at
scale. A 413 rejects the request before creating jobs, so split the body or
upload file content separately. `DELETE /delete-job` is a **soft delete** — it
hides the job and leaves stored result files intact.
