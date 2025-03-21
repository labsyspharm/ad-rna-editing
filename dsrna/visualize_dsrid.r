library(tidyverse)


x <- read_tsv("~/Downloads/dsRID_prediction_randomf.tsv")


x_short <- x %>%
  select(
    std_start, std_end,
    mean_start, mean_end,
    len_skip, coverage, num_skip,
    skip_ratio,
    chr, start, end,
    pred_1
  )
