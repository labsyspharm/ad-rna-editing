library(tidyverse)
library(synExtra)
library(data.table)

synapser::synLogin()
syn <- synDownloader("~/data", .cache = TRUE)

rm_tar_path <- syn("syn55256520")

untar(
  rm_tar_path,
  exdir = "repeatmasker_raw"
)

rm_raw_files <- list.files(
  "repeatmasker_raw/Human.RepeatMasker.map", full.names = TRUE, pattern = "chr.+out$"
) %>%
  set_names(
    str_match(., "(chr[0-9XYM]+)\\.fa\\.out")[, 2]
  )

rm_raw <- map(
  rm_raw_files,
  \(x) read_table(
    x,
    col_names = c(
      "sw_score", "perc_div", "perc_del", "perc_ins",
      "chromosome", "start", "end", "left", "strand",
      "rep_name", "class_family", "rep_start", "rep_end",
      "rep_left", "id", "dummy"
    ),
    skip = 2
  )
)

raw_concat <- bind_rows(rm_raw)

write_csv(
  raw_concat,
  "repeatmasker_raw/repeatmasker_raw.csv.gz"
)

synStoreMany(
  "repeatmasker_raw/repeatmasker_raw.csv.gz",
  parentId = "syn55256519",
  forceVersion = FALSE,
  used = "https://www-girinst-org.ezp-prod1.hul.harvard.edu/downloads/repeatmaskedgenomes/"
)
