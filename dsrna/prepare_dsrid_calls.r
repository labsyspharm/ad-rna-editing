library(tidyverse)
library(synExtra)
library(data.table)

synapser::synLogin()
syn <- synDownloader("~/data", .cache = TRUE)

dsrid_raw_files <- Sys.glob(
  "*/dsRID_prediction_randomf.tsv"
) %>%
  str_subset("pred_test", negate = TRUE) %>%
  tibble(file_path = .) %>%
  mutate(
    sample_id = str_match(
      file_path, "^(.+)\\.fastq\\.gz"
    )[, 2],
    chromosome = str_match(
      file_path, "(chr[0-9XYM]+)/"
    )[, 2],
    raw = map(
      file_path,
      read_tsv
    )
  )

dsrid_raw_concat <- dsrid_raw_files %>%
  group_by(sample_id) %>%
  summarize(
    raw = bind_rows(raw) %>%
      list(),
    .groups = "drop"
  ) %>%
  mutate(
    file_path = paste0(sample_id, "_dsRID_prediction_randomf.csv.gz")
  )

pwalk(
  dsrid_raw_concat,
  \(raw, file_path, ...) {
    write_csv(
      raw,
      file_path
    )
  }
)

synStoreMany(
  dsrid_raw_concat$file_path,
  parentId = "syn55256638",
  forceVersion = FALSE
)
