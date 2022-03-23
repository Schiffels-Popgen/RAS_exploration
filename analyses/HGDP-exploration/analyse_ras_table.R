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

ras_table %>% dplyr::filter(dataset == "1000G" & rasAF == "Common" & grepl("<?????A>", Left)) %>%
  dplyr::select(Left, Right, RAS, StdErr) %>%
  tidyr::pivot_wider(names_from = Right, values_from = c(RAS, StdErr)) %>%
  dplyr::mutate(
    ratio = RAS_FIN2 / RAS_IBS2,
    ratioErr = sqrt((StdErr_FIN2 / RAS_IBS2) ^ 2 + (StdErr_IBS2 * RAS_FIN2 / (RAS_IBS2^2))^2)
  ) %>% dplyr::select(Left, ratio, ratioErr) %>%
  ggplot(aes(x = Left, y = ratio)) + geom_point()
    # geom_errorbar(aes(ymin = ratio - ratioErr, ymax = ratio + ratioErr))

