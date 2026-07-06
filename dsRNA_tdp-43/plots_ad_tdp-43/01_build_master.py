"""
01_build_master.py
------------------
to build the specimen level master metadata table and the TDP-43 extract, and
write them to OUTDIR
"""
import os
import numpy as np
import pandas as pd
import data_prep as dp

OUTDIR = os.environ.get("ROSMAP_OUTDIR", "figures")
os.makedirs(OUTDIR, exist_ok=True)




def main():
    # --- TDP-43 extract here ---
    cs = pd.read_excel(dp._p("radc_cross"))[["projid", "study", "tdp_st4"]]
    tdp = cs.loc[cs["tdp_st4"].notna()].copy()
    tdp["tdp_st4"] = tdp["tdp_st4"].astype(int)
    tdp["tdp_st4_stage"] = tdp["tdp_st4"].map(
        {0: "None", 1: "Amygdala", 2: "Amygdala + Limbic",
         3: "Amygdala + Limbic + Neocortical"})
    tdp = tdp.sort_values("projid")
    tdp.to_csv(f"{OUTDIR}/tdp43_extract.csv", index=False)
    print("tdp43_extract:", tdp.shape)

    #### master file
    master = dp._derive_common(dp.build_master())
    keep = ["specimenID", "brain_region", "individualID", "projid", "tdp43_match",
            "n_expressed", "STMN2short", "UNC13A_CE1", "UNC13A_CE2",
            "class_low", "class_high", "tdp_st4", "tdp_st4_stage",
            "tdp_beyond_amygdala", "Study", "msex", "age_death", "apoe_genotype",
            "braaksc", "ceradsc", "cogdx", "pmi", "educ", "race", "RIN",
            "classification_file"]
    keep = [c for c in keep if c in master.columns]
    master[keep].to_csv(f"{OUTDIR}/ROSMAP_TDP43_master_metadata.csv", index=False)
    matched = master[master.tdp43_match == "matched"]
    matched[keep].to_csv(f"{OUTDIR}/ROSMAP_TDP43_matched_subset.csv", index=False)
    print("master:", master.shape, "| matched:", matched.shape)

    # summary
    for r in ["HCN", "PCC", "DLPFC"]:
        g = master[master.brain_region == r]
        print(f"  {r}: n={len(g)}  linked={g.projid.notna().sum()} "
              f"matched={(g.tdp43_match == 'matched').sum()}")


if __name__ == "__main__":
    main()
