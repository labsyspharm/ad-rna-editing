# MAJIQ → crypTE-Exon pipeline (Liu et al., 14 sorted samples)

Runs **MAJIQ** on the Liu RNA-seq to detect **local splicing variations (LSVs)**,
then flags the junctions whose splice site falls inside a **transposable element**
— i.e. the **crypTE-Exon** class from Bolger et al. 2026 (Fig 4).

MAJIQ reads junctions **directly from BAMs**, including de novo/cryptic ones, so
**STAR `SJ.out.tab` files are not needed.**

## Scope — which crypTE classes this covers
| crypTE class | Definition (paper) | Covered here? |
|---|---|---|
| **Exon** | splice donor/acceptor **inside** a TE | **Yes** — this is exactly what LSV junctions give |
| TSS | TE ≤50 bp from transcript **start** | No — needs transcript 5′ ends |
| ApA | TE ≤50 bp from transcript **end** | No — needs transcript 3′ ends |

MAJIQ is junction-centric, so it delivers the **Exon** class robustly. TSS/ApA are
transcript-**end** classes and require full transcript models (StringTie/bambu);
they are deliberately out of scope for a junction-based tool.

## What you must supply (the two things I can't)
1. **MAJIQ academic license file** — accept the license at <https://majiq.biociphers.org/>
   and download it; set `MAJIQ_LICENSE_FILE`. MAJIQ will not run without it.
2. **The 14 Liu sorted BAMs, indexed**, named by Run ID:
   `SRR8571937.bam` + `SRR8571937.bam.bai`, … (see list below). If yours are
   named/sorted differently:
   ```bash
   samtools sort -@8 -o SRR8571937.bam your_input.bam && samtools index SRR8571937.bam
   ```
   Put them all in one directory and set `BAMDIR`.

## The 14 sorted samples (2 unsorted dropped, per your instruction)
`liu_sample_groups.csv` — grouped for the deltapsi contrast **TDPneg − TDPpos**
(positive dPSI = higher inclusion on TDP-43 loss, the cryptic direction):

| group | n | Run IDs |
|---|---|---|
| **TDPneg** (grp1) | 7 | SRR8571938, 940, 942, 945, 948, 950, 952 |
| **TDPpos** (grp2) | 7 | SRR8571937, 939, 941, 944, 947, 949, 951 |

## Run it
```bash
bash 00_install_majiq.sh            # one-time: conda env 'majiq' + htslib + MAJIQ
conda activate majiq
export MAJIQ_LICENSE_FILE=/path/to/your_academic.lic
# edit params.sh: BAMDIR, REF_GFF3 (or let 01 fetch GENCODE v47), RMSK_CSV, STRANDNESS
bash run_all.sh                     # build -> deltapsi -> voila tsv -> crypTE
```
Or step by step: `01_build.sh` → `02_deltapsi.sh` → `03_voila_tsv.sh` → `04_crypte_from_lsv.py`.

Point `RMSK_CSV` at this project's **`repeatmasker_raw.csv`** (or any RepeatMasker
table with the same columns; a 6-col BED also works via `--rmsk-bed`).

## Two settings that silently break results if wrong
- **Strandedness.** Liu used the NuGEN **Ovation RNA-Seq System V2 → non-directional**,
  so `STRANDNESS="None"`. If unsure, run MAJIQ's builder on one sample each way and
  keep the setting that yields non-trivial junction coverage. A wrong value zeroes coverage.
- **Contig naming.** The GFF3 and the BAMs must both be `chr1…` or both `1…`.
  `01_build.sh` checks this and aborts on a mismatch. (`04` normalizes the `chr`
  prefix itself, and handles RepeatMasker's `C` = minus strand.)

## Outputs (`results/crypte/`)
- **`crypte_exon_candidates.tsv`** — de novo junctions, `P(changing) ≥ 0.90`,
  `|dPSI| ≥ 0.10`, up in TDPneg, with a splice site inside a TE. The headline result.
- **`crypte_exon_junctions.tsv.gz`** — every junction×TE hit (unfiltered), with
  `de_novo`, `dPSI`, `prob_changing`, TE name/class/family/divergence, and TE
  **orientation** vs the host strand (antisense 3′SS = classic exonization; sense =
  possible polyA/readthrough — the paper does not separate these).
- **`crypte_summary.txt`** — counts by TE class, orientation, and top families.

## Caveats
- `deltapsi` compares two groups of 7; it does **not** model the 7-patient pairing.
  You can exploit the pairing downstream (e.g. require a within-patient PSI shift).
- These are **candidate** crypTE-Exons: MAJIQ + TE-overlap. Bolger additionally
  require the TE portion not to overlap an annotated exon — add that filter against
  your GFF3 if you want to match their stringency exactly.
- Pipeline is validated on synthetic + real-schema fixtures; it has **not** been run
  on the Liu BAMs here because those BAMs and the license are not in this workspace.
