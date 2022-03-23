library(magrittr)
library(ggplot2)

read_ras_table <- function(full_filename) {
  fn <- basename(full_filename)
  m <- stringr::str_match(fn, "AncientBritish_(HGDP|1000G)_ras([0-9A-Za-z]+).table.txt")
  dataset <- m[2]
  af <- m[3]
  readr::read_tsv(full_filename, col_types="ccidd") %>%
    dplyr::mutate(dataset = dataset, rasAF = af)
}

ras_table <- list.files("analyses/HGDP-exploration", pattern = "^AncientBritish", full.names=TRUE) %>%
  purrr::map_dfr(~read_ras_table(.))
  
  readr::read_tsv("analyses/HGDP-exploration/AncientBritish_HGDP_ras.table.txt")

