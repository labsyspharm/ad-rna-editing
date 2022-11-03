library(tidyverse)
library(synExtra)
library(synapser)

# Select MAYO samples for download. here we don't really care so much
# about rRNA depletion vs poly-A enrichment because the main reason for
# looking at MAYO data is that we have gold standard TDP-43 pathology
# annotation for validating our STMN2 splicing predictor.

synLogin()

syn <- synDownloader("~/data", .cache = TRUE)

clinical_meta <- syn("syn14031984") %>%
  read_tsv()

fastq_meta <- bind_rows(
  synChildren("syn8612213") %>%
    enframe("file", "syn_id") %>%
    bind_cols(
      str_match(
        .$file,
        r"{(.*)\..*\.r[12]\.fastq\.gz}"
      ) %>%
        magrittr::set_colnames(
          c("file", "specimen_id")
        ) %>%
        as_tibble() %>%
        select(-file)
    ),
  synChildren("syn8612203") %>%
    enframe("file", "syn_id") %>%
    bind_cols(
      str_match(
        .$file,
        r"{(.*)\.r[12]\.fastq\.gz}"
      ) %>%
        magrittr::set_colnames(
          c("file", "specimen_id")
        ) %>%
        as_tibble() %>%
        select(-file)
    )
)

set.seed(42)
fastq_subset <- fastq_meta %>%
  inner_join(
    clinical_meta %>%
      drop_na(),
      # drop_na(`TDP-43`) %>%
      # group_by(`TDP-43`) %>%
      # slice_sample(n = 30) %>%
      # ungroup(),
    by = c("specimen_id" = "NETdbID")
  )

write_csv(
  fastq_subset, "selected_mayo_samples.csv"
)

