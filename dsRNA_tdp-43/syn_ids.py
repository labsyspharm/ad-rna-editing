#!/usr/bin/env python3
"""
syn_ids.py -- Synapse IDs for the Liu crypTE plots data (single source of truth).

All synIDs used by the plotting scripts live here so the mapping is reviewable in
one place (mirrors how the lab's R scripts hard-code syn("synXXXX") at the top of
each chunk). Folder layout on Synapse:

    inputs/        syn78078325
    intermediates/ syn78078328
    outputs/       syn78078327
"""

# ---- parent folders (upload targets for store()) ----------------------------
FOLDERS = {
    "inputs":        "syn78078325",
    "intermediates": "syn78078328",
    "outputs":       "syn78078327",
}

# ---- input files (downloaded by ID) -----------------------------------------
INPUTS = {
    "cryptic_junctions_TE_annotated_14samples.csv.gz": "syn78078345",
    "de_te_locus_14.csv.gz":                           "syn78078393",
    "full_AS_pattern_classification_14samples.csv":    "syn78078337",
    "liu_te_annotation.csv.gz":                        "syn78078384",
    "retained_introns_TE_annotated_14samples.csv.gz":  "syn78078344",
}

# ---- intermediate (output of script 01, input of script 02) -----------------
INTERMEDIATES = {
    "cryptic_junctions_AS_pattern_14samples.csv.gz":   "syn78078381",
}

# ---- existing output entities (re-upload creates a new VERSION of these) -----
# Figures (Fig_excess_TE_misplicing_refined_14samples.png,
# Fig_AS_taxonomy_modulizer_names.png) are created new in outputs/ on first run.
OUTPUTS = {
    "excess_TE_misplicing_refined_events_14samples.csv": "syn78078338",
}
