"""
data_prep.py
------------
Load the ROSMAP + RADC source files, then build the specimen-level master metadata
table, and also derive every analysis variable used across the figure scripts

    ROSMAP_clinical__2_.csv
    ROSMAP_biospecimen_metadata__1_.csv
    ROSMAP_assay_rnaSeq_metadata__1_.csv
    rosmap_classification__1_.csv
    rosmap_metadata__3_.csv
    dataset_1698_cross-sectional_06-10-2026__1_.xlsx   (info for autopsy/clinical)
"""
import os
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

DATA_DIR = os.environ.get("ROSMAP_DATA_DIR", "/data/uploads")

FILES = {
    "clinical":       "ROSMAP_clinical__2_.csv",
    "biospecimen":    "ROSMAP_biospecimen_metadata__1_.csv",
    "assay":          "ROSMAP_assay_rnaSeq_metadata__1_.csv",
    "classification": "rosmap_classification__1_.csv",
    "metadata":       "rosmap_metadata__3_.csv",
    "radc_cross":     "dataset_1698_cross-sectional_06-10-2026__1_.xlsx",
}

# TDP-43 stage labels (tdp_st4)
TDP_STAGE = {0: "None (no regions)", 1: "Amygdala",
             2: "Amygdala + Limbic", 3: "Amygdala + Limbic + Neocortical"}
DX_MAP = {1: "NCI", 2: "MCI", 3: "MCI", 4: "AD", 5: "AD", 6: "Other"}


def _p(key):
    return os.path.join(DATA_DIR, FILES[key])


def toint(s):
    return pd.to_numeric(s, errors="coerce").astype("Int64")


def load_sources():
    clin = pd.read_csv(_p("clinical"), dtype=str)
    bios = pd.read_csv(_p("biospecimen"), dtype=str)
    assay = pd.read_csv(_p("assay"), dtype=str)
    cls = pd.read_csv(_p("classification"), dtype=str)
    meta = pd.read_csv(_p("metadata"), dtype=str)
    return clin, bios, assay, cls, meta


def build_master():
  
    clin, bios, assay, cls, meta = load_sources()
    radc = pd.read_excel(_p("radc_cross"))[["projid", "tdp_st4"]]
    for d in (clin, meta):
        d["projid_int"] = toint(d["projid"])
    radc["projid_int"] = toint(radc["projid"])
    radc_tdp = radc.dropna(subset=["tdp_st4"]).drop_duplicates("projid_int")

    meta_map = meta.drop_duplicates("sample_id").set_index("sample_id")
    bios_map = (bios.dropna(subset=["individualID"])
                    .drop_duplicates("specimenID")
                    .set_index("specimenID")["individualID"])
    ind2proj = (clin.dropna(subset=["individualID"])
                    .drop_duplicates("individualID")
                    .set_index("individualID")["projid_int"])

    df = cls.rename(columns={"UNC13A-CE1": "UNC13A_CE1",
                             "UNC13A-CE2": "UNC13A_CE2",
                             "file": "classification_file"}).copy()

    # resolve individualID + projid
    df["individualID"] = (df["specimenID"].map(meta_map["individualID"])
                          .fillna(df["specimenID"].map(bios_map)))
    df["projid_int"] = (df["specimenID"].map(meta_map["projid_int"])
                        .fillna(df["individualID"].map(ind2proj)))

    # autopsy TDP-43 stage
    df = df.merge(radc_tdp[["projid_int", "tdp_st4"]], on="projid_int", how="left")
    df["tdp_st4"] = pd.to_numeric(df["tdp_st4"], errors="coerce").astype("Int64")

    # clinical / neuropath by projid
    cl = clin.drop_duplicates("projid_int").set_index("projid_int")
    for c in ["Study", "msex", "educ", "race", "apoe_genotype", "age_death",
              "braaksc", "ceradsc", "cogdx", "pmi"]:
        df[c] = df["projid_int"].map(cl[c])

    am = assay.drop_duplicates("specimenID").set_index("specimenID")
    co = lambda p, f: p.where(p.notna() & (p.astype(str) != "nan"), f)
    df["RIN"] = co(df["specimenID"].map(meta_map["rin"]),
                   df["specimenID"].map(am["RIN"]))

    df["tdp43_match"] = np.where(df["tdp_st4"].notna(), "matched", "unmatched")
    df["projid"] = df["projid_int"].astype("string")
    return df


def _derive_common(df):
    
    df = df.copy()
    for c in ["tdp_st4", "n_expressed", "braaksc", "ceradsc", "cogdx",
              "age_death", "msex", "RIN", "pmi"]:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # cryptic-exon booleans 
    for c in ["STMN2short", "UNC13A_CE1", "UNC13A_CE2"]:
        if c in df:
            df[c + "_b"] = df[c].astype(str).str.upper().map({"TRUE": 1, "FALSE": 0})
    if "n_expressed" in df:
        df["molpos"] = (df["n_expressed"] >= 1).astype("float")

    ### autopsy TDP-43 derivations
    if "tdp_st4" in df:
        df["tdp_st4_stage"] = df["tdp_st4"].map(TDP_STAGE)
        df["tdp_beyond_amygdala"] = df["tdp_st4"].map({0: 0, 1: 0, 2: 1, 3: 1})

    # clinical diagnosis groups
    if "cogdx" in df:
        df["dx"] = df["cogdx"].map(DX_MAP)

    
    ## Braak groupings
    def g3(v):
        if pd.isna(v):
            return "NA"
        v = int(v)
        return "Braak 0–II" if v <= 2 else ("Braak III–IV" if v <= 4 else "Braak V–VI")

    def g4(v):  
        if pd.isna(v):
            return "NA"
        v = int(v)
        return "0–III" if v <= 3 else {4: "IV", 5: "V", 6: "VI"}[v]

    def gbin(v):
        if pd.isna(v):
            return np.nan
        return "Control (Braak 0–III)" if v <= 3 else "Case (Braak IV–VI)"

    if "braaksc" in df:
        df["braak_grp3"] = df["braaksc"].map(g3)
        df["braak_grp4"] = df["braaksc"].map(g4)
        df["braak_bin"] = df["braaksc"].apply(gbin)

    def a4(v):
        if pd.isna(v):
            return np.nan
        try:
            return 1.0 if "4" in str(int(float(v))) else 0.0
        except Exception:
            return np.nan
    if "apoe_genotype" in df:
        df["apoe4b"] = df["apoe_genotype"].apply(a4)
    return df


def pcc_master():
    """PCC specimens only (n=546 and one row per participant)."""
    m = _derive_common(build_master())
    return m[m.brain_region == "PCC"].copy()


def radc_full():
    
    cs = _derive_common(pd.read_excel(_p("radc_cross")))
    cs["AD"] = (cs["cogdx"].isin([4, 5])).astype(float)
    cs.loc[cs.cogdx.isna(), "AD"] = np.nan
    cs["beyond"] = (cs["tdp_st4"] >= 2).astype(float)
    cs.loc[cs.tdp_st4.isna(), "beyond"] = np.nan
    cs["male"] = cs["msex"]
    cs["cerad_plaque"] = 4 - cs["ceradsc"]           
    cs["age10"] = cs["age_death"] / 10.0             
    cs["HS"] = pd.to_numeric(cs.get("hip_scl_yn_mid"), errors="coerce")
    return cs


if __name__ == "__main__":
    m = pcc_master()
    print("PCC master:", m.shape)
    print(m[["specimenID", "projid", "tdp_st4", "n_expressed",
             "braak_grp3", "dx"]].head())
