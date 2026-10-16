# MAJIQ → crypTE pipeline — Synapse edition

The MAJIQ → crypTE-Exon pipeline, wired to pull inputs from Synapse and push
outputs back with provenance. The MAJIQ/crypTE **compute is unchanged** — Synapse
I/O is added in the bash orchestration layer via the `synapse` CLI (the same
convention as `job_scripts/download_synapse.sh`), so the pure-compute scripts stay
portable.

## What changed vs the local pipeline
| File | Change |
|---|---|
| `params.sh` | + Synapse block: `USE_SYNAPSE` toggle and all synIDs |
| `synapse_io.sh` | **new** — `syn_get` / `syn_store` wrappers over the `synapse` CLI |
| `00_install_majiq.sh` | + `pip install synapseclient` (provides the `synapse` CLI) |
| `01_build.sh` | fetches `liu_sample_groups.csv`; after build, stores the `.sj` bundle |
| `03_voila_tsv.sh` | stores the LSV deltapsi tsv (under its Synapse name) |
| `run_all.sh` | fetches `repeatmasker_raw.csv.gz` before step 04; stores the 3 crypTE outputs |
| `make_config.py`, `02_deltapsi.sh`, `04_crypte_from_lsv.py` | **byte-identical** (pure compute; their Synapse I/O is handled by the orchestration) |

## Authentication (once)
The `synapse` CLI reads `SYNAPSE_AUTH_TOKEN` or `~/.synapseConfig`. Run
`synapse config` once, or `export SYNAPSE_AUTH_TOKEN=...`. No token is in the code.
Set `USE_SYNAPSE=0` in `params.sh` (or the environment) to run fully locally.

## Synapse layout (project folder `MAJIQ / majiq_liu_pipeline_data`, syn78074953)
```
inputs/        syn78074955   liu_sample_groups.csv (syn78074967), repeatmasker_raw.csv.gz (syn78075087)
intermediates/
  build/       syn78074959   liu_sj_files.tar.gz (syn78075826), sj_SRR8571952_14sample.tar.gz (syn78074999)
  voila/       syn78074961   TDPneg_TDPpos_LSV_deltapsi.tsv (syn78075015)
outputs/       syn78074956   crypte_exon_candidates_annotated.csv (syn78074966),
                             crypte_exon_junctions.tsv.gz (syn78074972), crypte_summary.txt (syn78074969)
```

## Synapse I/O per step
| Step | fetches | stores (→ folder) |
|---|---|---|
| `01_build.sh` | `liu_sample_groups.csv` (syn78074967) | `liu_sj_files.tar.gz` → build/ (syn78074959) |
| `02_deltapsi.sh` | — (uses local `.majiq` from 01) | — (the `.deltapsi.voila` binary stays local) |
| `03_voila_tsv.sh` | — (local splicegraph + `.voila`) | `…_LSV_deltapsi.tsv` → voila/ (syn78074961) |
| `04` (via `run_all.sh`) | `repeatmasker_raw.csv.gz` (syn78075087) | `crypte_exon_junctions.tsv.gz`, `crypte_exon_candidates.tsv`, `crypte_summary.txt` → outputs/ (syn78074956) |

Uploads target the parent **folder** ID, so re-running creates a new **version** of
a file that already exists there and a new entity otherwise. `used=` provenance
(the upstream synID) is attached to each upload.

NOT on Synapse (external / local): the raw BAMs, the GFF3 annotation, the MAJIQ
license, and the large MAJIQ binaries (`.majiq`, `splicegraph.sql`, `.deltapsi.voila`).
The `.sj` bundle on Synapse (syn78075826 / syn78074999) lets you resume the build
without the BAMs via `majiq build --incremental`.

## Run
```bash
bash 00_install_majiq.sh          # once (installs MAJIQ + synapseclient)
conda activate majiq
export SYNAPSE_AUTH_TOKEN=...      # or: synapse config
# edit params.sh: BAMDIR, MAJIQ_LICENSE_FILE, REF_GFF3 (USE_SYNAPSE=1 by default)
bash run_all.sh                    # build → deltapsi → voila tsv → crypTE, with Synapse I/O
```
