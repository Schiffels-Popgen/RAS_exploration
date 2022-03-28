library(magrittr)
library(ggplot2)

pat <- "AncientBritish_(HGDP|1000G)_1240K_ras([0-9A-Za-z]+).table.txt"

popNames <- tibble::tribble(
  ~Left, ~Group,
  "<I0156.SG>", "AncientBritish",
  "<I0157.SG>", "AncientBritish",
  "<I0159.SG>", "AncientBritish",
  "<I0160.SG>", "AncientBritish",
  "<I0161.SG>", "AncientBritish",
  "<I0769.SG>", "AncientBritish",
  "<I0773.SG>", "AncientBritish",
  "<I0774.SG>", "AncientBritish",
  "<I0777.SG>", "AncientBritish",
  "<I0789.SG>", "AncientBritish",
  "<3DT16.SG>",  "AncientBritish",
  "<3DT26.SG>",  "AncientBritish",
  "<6DT18.SG>",  "AncientBritish",
  "<6DT21.SG>",  "AncientBritish",
  "<6DT22.SG>",  "AncientBritish",
  "<6DT23.SG>",  "AncientBritish",
  "<6DT3.SG>",   "AncientBritish",
  "<M1489.SG>",  "AncientBritish",    
  "<NO3423.SG>", "AncientBritish",
  "<NA12889.SG>", "CEU",
  "<NA12890.SG>", "CEU",
  "<HG00383.SG>", "FIN",
  "<HG00384.SG>", "FIN",
  "<HG01791.SG>", "GBR",
  "<HG02215.SG>", "GBR",
  "<HG02238.SG>", "IBS", 
  "<HG02239.SG>", "IBS",
  "<NA20828.SG>", "TSI",
  "<NA20832.SG>", "TSI",
  "<NA18749.SG>", "CHB",
  "<NA18757.SG>", "CHB",
  "<NA19248.SG>", "YRI",
  "<NA19256.SG>", "YRI",
  "<NA19257.SG>", "YRI",
  "<HG02348.SG>", "PEL",
  "<HG02425.SG>", "PEL"
)

read_ras_table <- function(full_filename) {
  fn <- basename(full_filename)
  m <- stringr::str_match(fn, pat)
  dataset <- m[2]
  af <- m[3]
  readr::read_tsv(full_filename, col_types="ccidd") %>%
    dplyr::mutate(dataset = dataset, rasAF = af) %>%
    dplyr::left_join(popNames)
}

ras_table <- list.files("analyses/HGDP-exploration", pattern = pat, full.names=TRUE) %>%
  purrr::map_dfr(~read_ras_table(.))

ras_table %>% dplyr::filter(dataset == "1000G" & rasAF == "05" & !(Right %in% c("CHB2", "YRI2")) &
                              !(Group %in% c("YRI", "CHB", "PEL"))) %>%
  dplyr::select(Left, Group, Right, RAS, StdErr) %>%
  ggplot(aes(x = Left, y = RAS, col = Group)) + geom_point() +
  geom_errorbar(aes(ymin = RAS - StdErr, ymax = RAS + StdErr)) +
  facet_wrap(~Right) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

for(freq in c("01", "02", "05", "10", "20", "All", "Common")) {
  plot <- ras_table %>% dplyr::filter(dataset == "1000G" & rasAF == freq & !(Right %in% c("CHB2", "YRI2")) &
                                        !(Group %in% c("YRI", "CHB", "PEL"))) %>%
    dplyr::select(Left, Group, Right, RAS, StdErr) %>%
    ggplot(aes(x = Left, y = RAS, col = Group)) + geom_point() +
    geom_errorbar(aes(ymin = RAS - StdErr, ymax = RAS + StdErr)) +
    facet_wrap(~Right) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  fn <- paste0("analyses/HGDP-exploration/scatter_plot_1000G_1240K_ras", freq, ".pdf")
  ggsave(fn, plot)
}

# Checking coverage dependency
# Checking dependency on Coverage
ras_table %>%
  dplyr::filter(dataset == "1000G" & rasAF == "Common") %>%
  ggplot() + geom_point(aes(x = Norm, y = RAS, col = Right))
