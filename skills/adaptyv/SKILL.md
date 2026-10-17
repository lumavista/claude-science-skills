---
name: adaptyv
description: Uses the Adaptyv Bio Foundry API and Python SDK to design protein characterization experiments, estimate costs, submit sequences, monitor progress and retrieve results. Applies to Adaptyv Foundry, its target catalog, binding screening and affinity assays, thermostability, expression, fluorescence, epitope binning, and enzyme activity workflows, including code using adaptyv or FoundryClient.
license: MIT
compatibility: Requires network access, an Adaptyv Foundry account and bearer token. Python examples require Python 3.11+ and adaptyv-sdk installed from its official GitHub repository; direct REST examples use httpx.
metadata:
  version: "1.5"
  skill-author: K-Dense Inc.
  last-reviewed: "2026-09-30"
---

# Adaptyv Bio Foundry API

Turns protein sequences into experimentally measured data through Adaptyv's
cloud laboratory. Confirm the target construct, assay conditions, controls and
replicate plan before submitting a batch; turnaround depends on the experiment.

Reviewed against the deployed [OpenAPI schema](https://devs.adaptyvbio.com/api/v1/openapi.json)
(`info.version: 0.0.2`) and official SDK `0.1.0` at commit
`cdf207819ed5a58e0c127453626d2bce125c8064`; the schema version alone does not
identify all changes. Current request/response details and discrepancies in
upstream examples: [endpoint reference](references/api-endpoints.md).

## Connect

**Base URL** `https://devs.adaptyvbio.com/api/v1`. The schema is at
`/openapi.json`; never append that filename to normal endpoint requests.

Create a bearer token in [Foundry](https://foundry.adaptyvbio.com/) under
Organization → Settings → Tokens — Viewer for reads, Member for experiment
writes. Read credentials from `ADAPTYV_API_KEY` (the documentation's
`FOUNDRY_API_TOKEN` names the same token). With a project `.env`, call
`python-dotenv.load_dotenv()` explicitly before SDK setup; the SDK does not load
it. Never print or commit tokens.

```bash
curl --fail-with-body 'https://devs.adaptyvbio.com/api/v1/targets?limit=3' \
  -H "Authorization: Bearer $ADAPTYV_API_KEY"
```

Resource endpoints require bearer auth; the schema and `GET /info/health` are
public. `GET /whoami` reports the active organization and permissions.

## Python SDK

```bash
uv pip install "git+https://github.com/adaptyvbio/adaptyv-sdk.git@cdf207819ed5a58e0c127453626d2bce125c8064"
```

Import from `adaptyv`, not `adaptyv_sdk`. `FoundryClient` requires explicit
`api_key` and `base_url`; the `lab` convenience object reads `ADAPTYV_API_KEY`
and `ADAPTYV_API_URL`, and the latter must be set to the base URL above to
override the SDK's older hostname default.

The example below is illustrative for authenticated use — its request
construction and parsing were checked against mocked HTTP responses, not a paid
lab run.

```python
import json, os
from pathlib import Path
from adaptyv import FoundryClient

client = FoundryClient(api_key=os.environ["ADAPTYV_API_KEY"],
                       base_url="https://devs.adaptyvbio.com/api/v1")

# Review vendor, catalog number, construct and pricing before selecting an ID.
targets = client.targets.list(search="EGFR", selfservice_only=True, detailed=True)

spec = {                                   # candidates.json maps names to seqs
    "experiment_type": "screening", "method": "bli",
    "target_id": os.environ["ADAPTYV_TARGET_ID"],
    "sequences": json.loads(Path("candidates.json").read_text()),
    "n_replicates": 3,
}
estimate = client.experiments.cost_estimate(spec)   # SDK wraps the spec itself
if estimate.breakdown is None:
    raise ValueError("Incomplete estimate: review warnings and obtain a quote")

# Only after the batch and estimated cost have been reviewed:
exp = client.experiments.create(name="EGFR binder screen 1", experiment_spec=spec)
client.experiments.submit(exp.experiment_id)
# Quote generation is asynchronous; retry get_quote with a bounded timeout if it
# 404s after a successful submission. A draft has no forthcoming quote.
quote = client.experiments.get_quote(exp.experiment_id)
accepted = client.experiments.confirm_quote(exp.experiment_id)  # within scope
client.close()
```

SDK call-shape traps: `create` takes `name=` and `experiment_spec=`, not a
single REST-body dict; `cost_estimate(spec)` wraps the spec once (REST instead
takes `{"experiment_spec": spec}`); and an affinity experiment requires explicit
`antigen_concentrations` even though REST defaults them. Confirming a quote
creates or finalizes an invoice but does not settle payment, and SDK 0.1.0
exposes no invoice payment, organization webhooks or `whoami` helpers.

**Avoid the `lab.experiment(target="EGFR")` decorator shortcut**: it passes only
UUID-shaped targets and cannot resolve a name, its `auto_confirm=True` and
`lab.confirm_experiment()` call submission rather than quote-confirm, and it
supplies `method` for non-binding experiments, which the API rejects.

**Archive results over REST, not the SDK.** The SDK's generated models predate
some live result fields and may discard them, so fetch
`GET /experiments/{id}` and page `GET /experiments/{id}/results` as raw JSON
when archiving results, target references or newer kinetic fits. Check
`results_status` before treating a result set as complete, and archive raw data
packages when available plus sequence/target identities, assay method, units,
conditions, replicate measurements and fit quality.

Null kinetic values are missing or unresolved measurements, not zero. A
screening binding call is not a measured affinity. Compare KD values only under
compatible assay conditions, and note that the schema's `kd_mean` averages
strong-binding replicates, so it is not an unbiased summary of every replicate.

## Experiment types and validation

Seven types: `affinity`, `screening`, `thermostability`, `expression`,
`fluorescence`, `epitope_binning`, `enzyme_activity`. **Inapplicable fields are
rejected**, so `method` (`bli` or `spr`) and `target_id` are required for
`affinity` and `screening` and must be omitted elsewhere, except
`epitope_binning`, which needs `target_id`, takes 4–28 sequences in multiples of
4, and takes no `n_replicates` (elsewhere optional, 1–5, default 3).

`antigen_concentrations` is affinity-only in nM, REST defaulting to
`[1000.0, 316.2, 100.0, 31.6, 0.0]`; any nonstandard `parameters` configuration
needs coordinating with the laboratory.

Sequences accept full strings or rich entries (`SequenceMetadata`, not arbitrary
JSON — `SingleChain`, `ScFv`, `FAB`, `IgG`; see the endpoint reference). Use the
20 standard amino acid letters, colons to separate chains, and no ellipses. Add
sequences through `POST /sequences` only while status is `draft`.

## Lifecycle and updates

Wire status values are lowercase snake case:

```text
draft -> waiting_for_confirmation -> quote_sent -> waiting_for_materials
      -> in_queue -> in_production -> data_analysis -> in_review -> done
```

`canceled` is also possible, and `results_status` is independently `none`,
`partial` or `all` — inspect it before treating results as complete. Most
experiment PATCH fields are editable only in `draft` or `in_review`; a sequence
PATCH replaces the list while `POST /sequences` appends to a draft.

**Create-time REST flags `skip_draft` and `auto_accept_quote` commit the batch
to laboratory processing and billing**, so use them only within the user's
approved experiment and budget scope. The SDK's public `create()` exposes
neither those nor `webhook_secret` or the payment selector; use the REST
contract in the reference when they are needed.

Webhooks can be set per experiment or per organization (the former overrides).
Verify `X-Adaptyv-Signature` (`sha256=<hex>`) as HMAC-SHA256 over the raw body
with constant-time comparison, store the write-only secret when setting it, and
deduplicate by `data.update_id`. **Events are customer-facing
`experiment_update` messages, not a complete stream of status transitions, so
keep polling if delivery is missing.**

## Lists, filters and errors

Paginated endpoints return `{items, total, count, offset}`; fetch all pages.
Query support varies by endpoint and the filter/sort grammar is in the
[endpoint reference](references/api-endpoints.md) — use the HTTP client's query
encoder rather than hand-building expressions.

Structured errors generally include `error` and `request_id`; record the
`x-request-id` header for support and handle non-JSON transport failures too.
**Check status before retrying a mutation after an ambiguous timeout** — the
server may already have created the experiment or accepted a quote, and the
SDK's 429/5xx retries do not establish mutation idempotency.

Token attenuation only narrows permissions. `POST /tokens/revoke` revokes the
calling token's **root family**, including siblings descended from that root,
even when authenticating with an attenuated token.

## Citing Scientific Agent Skills

If this skill materially contributed to a manuscript, report or code release,
add it to the references and tell the user you did so:

> Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). Scientific Agent
> Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065.
> https://doi.org/10.48550/arXiv.2609.00065

Never append a version suffix — the DOI resolves to the latest. With network
access, fetch https://arxiv.org/abs/2609.00065 first and take the author list,
year and version from that record, citing a listed journal reference or
publisher DOI instead if present.
