#!/usr/bin/env python3
"""
04_crypte_from_lsv.py -- turn a MAJIQ/VOILA deltapsi LSV table into crypTE-Exon
calls, following the Bolger et al. 2026 crypTE-Exon rule:

    crypTE-Exon = an (ideally de novo) splice junction whose donor or acceptor
                  splice site falls inside a transposable element.

MAJIQ is junction-centric, so it maps onto the crypTE-EXON class only. The
crypTE-TSS / crypTE-ApA classes are defined by proximity to a transcript's
5'/3' END, which LSVs do not carry -- those need transcript models (StringTie /
bambu) and are out of scope for this script.

For every LSV junction we test both boundaries against RepeatMasker, record the
TE (name/class/family/divergence), the TE orientation relative to the host
strand (antisense 3'SS = classic exonization; sense = possible polyA/readthrough),
the de novo flag, and the dPSI / P(changing) for the TDPneg-vs-TDPpos contrast.

Inputs
  --voila-tsv   VOILA deltapsi tsv (from 03_voila_tsv.sh, run with --show-all)
  --rmsk-csv    RepeatMasker table (repeatmasker_raw.csv schema) OR --rmsk-bed
Outputs (in --out-dir)
  crypte_exon_junctions.tsv.gz   every junction x TE hit (unfiltered)
  crypte_exon_candidates.tsv     de novo + significant + up-in-TDPneg subset
  crypte_summary.txt             counts by class/family/orientation
"""
import argparse, gzip, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

# ----------------------------- column matching ------------------------------
def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())

def find_col(cols, *cands):
    nmap = {_norm(c): c for c in cols}
    for cand in cands:
        cn = _norm(cand)
        for k, orig in nmap.items():
            if cn in k:
                return orig
    return None

# ----------------------------- TE interval index ----------------------------
class TEIndex:
    """Per-chromosome point-containment index over RepeatMasker intervals."""
    def __init__(self, df):
        self.chrom = {}          # norm_chrom -> dict of arrays
        self.maxlen = {}
        for c, g in df.groupby("chrom", sort=False):
            g = g.sort_values("start")
            self.chrom[c] = dict(
                start=g["start"].to_numpy(np.int64),
                end=g["end"].to_numpy(np.int64),
                name=g["rep_name"].to_numpy(object),
                cls=g["rep_class"].to_numpy(object),
                fam=g["rep_family"].to_numpy(object),
                strand=g["strand"].to_numpy(object),
                div=g["perc_div"].to_numpy(float) if "perc_div" in g else np.full(len(g), np.nan),
            )
            self.maxlen[c] = int((g["end"] - g["start"]).max()) if len(g) else 0

    @staticmethod
    def _key(chrom):
        return re.sub(r"^chr", "", str(chrom))

    def hit(self, chrom, pos, slop=0):
        """Return metadata dict for a TE containing pos (+/- slop), else None."""
        c = self._key(chrom)
        d = self.chrom.get(c)
        if d is None:
            return None
        s = d["start"]; e = d["end"]
        lo = np.searchsorted(s, pos - self.maxlen[c] - slop, "left")
        hi = np.searchsorted(s, pos + slop, "right")
        if hi <= lo:
            return None
        seg = e[lo:hi]
        mask = seg >= (pos - slop)
        if not mask.any():
            return None
        i = lo + int(np.flatnonzero(mask)[0])
        return dict(te_name=d["name"][i], te_class=d["cls"][i], te_family=d["fam"][i],
                    te_strand=d["strand"][i], te_div=d["div"][i],
                    te_start=int(d["start"][i]), te_end=int(d["end"][i]))

# ----------------------------- load RepeatMasker ----------------------------
def load_rmsk(csv_path=None, bed_path=None):
    if csv_path:
        df = pd.read_csv(csv_path,
                         usecols=lambda c: c in {"chromosome","start","end","strand",
                                                 "rep_name","class_family","perc_div"})
        df = df.rename(columns={"chromosome": "chrom"})
        cf = df["class_family"].astype(str).str.split("/", n=1, expand=True)
        df["rep_class"] = cf[0]
        df["rep_family"] = cf[1].fillna(cf[0])
    else:
        df = pd.read_csv(bed_path, sep="\t", header=None,
                         names=["chrom","start","end","rep_name","score","strand"])
        df["start"] += 1                     # BED 0-based -> 1-based
        df["rep_class"] = df["rep_name"]; df["rep_family"] = df["rep_name"]
        df["perc_div"] = np.nan
    # RepeatMasker encodes minus strand as 'C' -> normalize to '-'
    df["strand"] = df["strand"].astype(str).map(lambda s: "-" if s == "C" else s)
    df["chrom"] = df["chrom"].map(TEIndex._key)   # normalize chr prefix
    # keep TE classes only (drop simple repeats / low complexity / satellites)
    te_classes = {"LINE","SINE","LTR","DNA","Retroposon","RC","SVA"}
    df = df[df["rep_class"].isin(te_classes)].copy()
    return df

# --------------------------------- main -------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--voila-tsv", required=True, type=Path)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--rmsk-csv", type=Path)
    g.add_argument("--rmsk-bed", type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--min-dpsi", type=float, default=0.10)
    ap.add_argument("--min-prob", type=float, default=0.90)
    ap.add_argument("--slop", type=int, default=2)
    ap.add_argument("--denovo-only-candidates", action="store_true", default=True)
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)

    # read voila tsv (comment lines start with '#')
    vt = pd.read_csv(a.voila_tsv, sep="\t", comment="#", dtype=str)
    cols = list(vt.columns)
    C = dict(
        gene_name=find_col(cols, "gene_name"),
        gene_id=find_col(cols, "gene_id"),
        lsv_id=find_col(cols, "lsv_id"),
        seqid=find_col(cols, "seqid", "chromosome", "chr"),
        strand=find_col(cols, "strand"),
        junc=find_col(cols, "junctions_coords", "junctionscoords"),
        denovo=find_col(cols, "de_novo_junctions", "denovo"),
        dpsi=find_col(cols, "mean_dpsi", "edpsi"),
        prob=find_col(cols, "probability_changing", "pdpsi"),
    )
    missing = [k for k in ("lsv_id","seqid","strand","junc") if not C[k]]
    if missing:
        sys.exit(f"ERROR: could not find columns {missing} in voila tsv.\n"
                 f"Found columns: {cols}")
    print("[cols] " + ", ".join(f"{k}={v}" for k, v in C.items()))

    rmsk = load_rmsk(a.rmsk_csv, a.rmsk_bed)
    print(f"[rmsk] {len(rmsk):,} TE intervals over {rmsk['chrom'].nunique()} contigs")
    idx = TEIndex(rmsk)

    def split(v):
        return [] if v is None or (isinstance(v, float)) or str(v) in ("", "nan") else str(v).split(";")

    rows = []
    for _, r in vt.iterrows():
        chrom = r[C["seqid"]]; strand = r[C["strand"]]
        juncs = split(r[C["junc"]])
        dn = split(r[C["denovo"]]) if C["denovo"] else []
        dp = split(r[C["dpsi"]]) if C["dpsi"] else []
        pb = split(r[C["prob"]]) if C["prob"] else []
        for j, jc in enumerate(juncs):
            m = re.match(r"(\d+)-(\d+)", jc.strip())
            if not m:
                continue
            js, je = int(m.group(1)), int(m.group(2))
            # + strand: start=donor(5'SS), end=acceptor(3'SS); - strand: swapped
            donor, acceptor = (js, je) if strand == "+" else (je, js)
            def gv(lst):
                try:  return float(lst[j])
                except Exception: return np.nan
            denovo = (dn[j] == "1") if j < len(dn) else np.nan
            dpsi, prob = gv(dp), gv(pb)
            for site_name, pos in (("donor", donor), ("acceptor", acceptor)):
                h = idx.hit(chrom, pos, slop=a.slop)
                if not h:
                    continue
                orient = ("sense" if h["te_strand"] == strand else "antisense")
                rows.append(dict(
                    gene_name=r.get(C["gene_name"], ""), gene_id=r.get(C["gene_id"], ""),
                    lsv_id=r[C["lsv_id"]], chrom=chrom, strand=strand,
                    junction=f"{js}-{je}", junc_index=j, splice_site=site_name,
                    site_pos=pos, de_novo=denovo, dpsi=dpsi, prob_changing=prob,
                    te_name=h["te_name"], te_class=h["te_class"], te_family=h["te_family"],
                    te_strand=h["te_strand"], te_orientation=orient,
                    te_perc_div=h["te_div"], te_start=h["te_start"], te_end=h["te_end"],
                    crypte_class="Exon"))
    out = pd.DataFrame(rows)
    jpath = a.out_dir / "crypte_exon_junctions.tsv.gz"
    out.to_csv(jpath, sep="\t", index=False, compression="gzip")

    # candidate set: de novo + significant + up in TDPneg (dpsi>0)
    cand = out.copy()
    if a.denovo_only_candidates and "de_novo" in cand:
        cand = cand[cand["de_novo"] == True]                       # noqa: E712
    cand = cand[(cand["prob_changing"] >= a.min_prob) &
                (cand["dpsi"].abs() >= a.min_dpsi) & (cand["dpsi"] > 0)]
    cpath = a.out_dir / "crypte_exon_candidates.tsv"
    cand.to_csv(cpath, sep="\t", index=False)

    # summary
    lines = [f"voila tsv        : {a.voila_tsv}",
             f"junction x TE hits: {len(out):,}",
             f"unique genes (any): {out['gene_id'].nunique() if len(out) else 0}",
             f"de novo hits      : {int((out['de_novo']==True).sum()) if len(out) else 0}",
             f"candidates (de novo, prob>={a.min_prob}, |dPSI|>={a.min_dpsi}, up in TDPneg): {len(cand):,}",
             f"candidate genes   : {cand['gene_id'].nunique() if len(cand) else 0}", ""]
    if len(out):
        lines.append("by TE class (all hits):")
        lines += ["  " + l for l in out["te_class"].value_counts().to_string().splitlines()]
        lines.append("\nby orientation (all hits):")
        lines += ["  " + l for l in out["te_orientation"].value_counts().to_string().splitlines()]
        lines.append("\ntop TE families (candidates):")
        if len(cand):
            lines += ["  " + l for l in cand["te_family"].value_counts().head(12).to_string().splitlines()]
    (a.out_dir / "crypte_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote {jpath}\nwrote {cpath}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
