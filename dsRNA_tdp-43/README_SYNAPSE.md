# Liu crypTE plots — Synapse edition

The three Liu crypTE plotting scripts, rewired to read inputs from Synapse by ID
and write outputs back to Synapse with provenance. The analysis/compute is
**unchanged** from the local versions — only file I/O is now Synapse-backed.

This mirrors the lab's R `synExtra` workflow (`synapser::synLogin()` →
`syn("synXXXX")` → `synStoreMany(..., parentId=, used=)`), implemented in Python
with `synapseclient`:

| R (`synExtra`) | Python (here, in `synapse_io.py`) |
|---|---|
| `synapser::synLogin()` | `synapse_io.login()` |
| `syn("synXXXX")` | `synapse_io.fetch("synXXXX")` → local path |
| `synStoreMany(p, parentId=, used=)` | `synapse_io.store(p, parent_id, used=, executed=)` |

## Authentication (do this once)
`synapseclient` reads your Synapse personal access token from either:
- the `SYNAPSE_AUTH_TOKEN` environment variable, or
- `~/.synapseConfig` (run `synapse config` once).

No token is stored in the code. Create a PAT at
https://www.synapse.org/ → Account Settings → Personal Access Tokens.

## Files
| File | Role |
|---|---|
| `syn_ids.py` | **all synIDs** (inputs / intermediates / output folders) — the single place to edit if IDs change |
| `synapse_io.py` | `login()` / `fetch()` / `store()` helpers over `synapseclient` |
| `01_crypte_AS_pattern_classify.py` | splicegraph + cryptic junctions → AS-pattern labels |
| `02_refined_excess_TE_missplicing.py` | refined excess-TE overlay (both splicing directions) |
| `03_AS_taxonomy_barplots.py` | Modulizer-name taxonomy bar plots |
| `requirements.txt` | deps (adds `synapseclient`) |

## Synapse I/O per script
| Script | fetches (inputs) | stores (→ folder) |
|---|---|---|
| **01** | `cryptic_junctions_TE_annotated_14samples.csv.gz` (syn78078345) | `cryptic_junctions_AS_pattern_14samples.csv.gz` → intermediates/ (syn78078328) |
| **02** | de (syn78078393), ann (syn78078384), C (syn78078345), IR (syn78078344), AS-pattern (syn78078381) | `Fig_excess_TE_misplicing_refined_14samples.png`, `excess_TE_misplicing_refined_events_14samples.csv` → outputs/ (syn78078327) |
| **03** | `full_AS_pattern_classification_14samples.csv` (syn78078337) | `Fig_AS_taxonomy_modulizer_names.png` → outputs/ (syn78078327) |

Uploads go to the parent **folder** ID, so re-running a script creates a new
**version** of an existing file (e.g. `excess_TE_..._events_14samples.csv`,
syn78078338) and a brand-new entity for a file that isn't there yet (the figures).
Provenance (`used=` input IDs, `executed=` script name) is attached to every upload.

## Run order
```bash
pip install -r requirements.txt
export SYNAPSE_AUTH_TOKEN=...          # or: synapse config

# 02 needs the AS-pattern intermediate (syn78078381), which already exists on
# Synapse — so 02 and 03 can be run directly:
python 02_refined_excess_TE_missplicing.py
python 03_AS_taxonomy_barplots.py

# 01 regenerates that intermediate. It additionally needs the MAJIQ splicegraph
# (results/build_off14/splicegraph.sql), which is a MAJIQ-pipeline output and is
# NOT in this Synapse folder — pass it as a local path:
python 01_crypte_AS_pattern_classify.py  /path/to/splicegraph.sql
```

Expected console output of 02:
`events overlaid: 305 (up 192, down 113) | expr-significant 11 (4%)`
`keep crypTE: up 74 down 44 | Alt-Intron: up 131 down 77`
