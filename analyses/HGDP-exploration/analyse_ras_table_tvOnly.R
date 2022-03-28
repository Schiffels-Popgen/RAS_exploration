library(magrittr)
library(ggplot2)

pat <- "AncientBritish_(HGDP|1000G)_ras([0-9A-Za-z]+)_TVonly.table.txt"

popNames <- tibble::tribble(
  ~Left, ~Group,
  "<12880A>", "AncientBritish",
  "<12881A>", "AncientBritish",
  "<12883A>", "AncientBritish",
  "<12884A>", "AncientBritish",
  "<12885A>", "AncientBritish",
  "<15558A>", "AncientBritish",
  "<15569A>", "AncientBritish",
  "<15570A>", "AncientBritish",
  "<15577A>", "AncientBritish",
  "<15579A>", "AncientBritish",
  "<3DT16>",  "AncientBritish",
  "<3DT26>",  "AncientBritish",
  "<6DT18>",  "AncientBritish",
  "<6DT21>",  "AncientBritish",
  "<6DT22>",  "AncientBritish",
  "<6DT23>",  "AncientBritish",
  "<6DT3>",   "AncientBritish",
  "<M1489>",  "AncientBritish",    
  "<NO3423>", "AncientBritish",
  "<NA12889>", "CEU",
  "<NA12890>", "CEU",
  "<HG00383>", "FIN",
  "<HG00384>", "FIN",
  "<HG01791>", "GBR",
  "<HG02215>", "GBR",
  "<HG02238>", "IBS", 
  "<HG02239>", "IBS",
  "<NA20828>", "TSI",
  "<NA20832>", "TSI",
  "<NA18749>", "CHB",
  "<NA18757>", "CHB",
  "<NA19248>", "YRI",
  "<NA19256>", "YRI",
  "<NA19257>", "YRI",
  "<HG02348>", "PEL",
  "<HG02425>", "PEL"
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

ras_table %>% dplyr::filter(dataset == "1000G" & rasAF == "01" & !(Right %in% c("CHB2", "YRI2")) &
                              !(Group %in% c("YRI", "CHB", "PEL"))) %>%
  dplyr::select(Left, Group, Right, RAS, StdErr) %>%
  ggplot(aes(x = Left, y = RAS, col = Group)) + geom_point() +
  geom_errorbar(aes(ymin = RAS - StdErr, ymax = RAS + StdErr)) +
  facet_wrap(~Right) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("sanalyses/HGDP-exploration/scatter_plot_1000G_ras01_tvOnly.pdf", plot1)

# No Fin plot
ras_table %>% dplyr::filter(dataset == "1000G" & rasAF == "01" & !(Right %in% c("CHB2", "YRI2")) &
                              !(Group %in% c("YRI", "CHB", "PEL", "FIN"))) %>%
  dplyr::select(Left, Group, Right, RAS, StdErr) %>%
  ggplot(aes(x = Left, y = RAS, col = Group)) + geom_point() +
  geom_errorbar(aes(ymin = RAS - StdErr, ymax = RAS + StdErr)) +
  facet_wrap(~Right) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Saving all noFin plots
for(freq in c("01", "02", "05", "10", "20", "All", "Common")) {
  plot <- ras_table %>% dplyr::filter(dataset == "1000G" & rasAF == freq & !(Right %in% c("CHB2", "YRI2")) &
                                !(Group %in% c("YRI", "CHB", "PEL", "FIN"))) %>%
    dplyr::select(Left, Group, Right, RAS, StdErr) %>%
    ggplot(aes(x = Left, y = RAS, col = Group)) + geom_point() +
    geom_errorbar(aes(ymin = RAS - StdErr, ymax = RAS + StdErr)) +
    facet_wrap(~Right) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  fn <- paste0("analyses/HGDP-exploration/scatter_plot_1000G_ras", freq, "_TVonly_noFin.pdf")
  ggsave(fn, plot)
}


# Checking F4-values
inds <- ras_table %>% dplyr::filter(dataset == "1000G") %>% dplyr::pull(Left) %>% unique()
rasAF <- ras_table %>% dplyr::filter(dataset == "1000G") %>% dplyr::pull(rasAF) %>% unique()
rasf4_table <- tidyr::expand_grid(Left1 = inds, Left2 = inds, rasAF = rasAF) %>%
  dplyr::mutate(dataset = "1000G") %>%
  dplyr::filter(Left1 != Left2) %>%
  dplyr::left_join(ras_table, by = c("Left1" = "Left", "dataset" = "dataset", "rasAF" = "rasAF")) %>%
  dplyr::rename(Norm1 = Norm, RAS1 = RAS, StdErr1 = StdErr, Group1 = Group) %>%
  dplyr::left_join(ras_table, by = c("Left2" = "Left", "dataset" = "dataset", "rasAF" = "rasAF", "Right" = "Right")) %>%
  dplyr::rename(Norm2 = Norm, RAS2 = RAS, StdErr2 = StdErr, Group2 = Group) %>%
  dplyr::mutate(
    f4name = paste0("RAS-F4(", Left1, ", ", Left2, ", ", Right, ", Outgroup)"),
    rasF4 = RAS1 - RAS2,
    rasF4err = sqrt(StdErr1^2 + StdErr2^2),
    rasF4Z = rasF4 / rasF4err
  )

rasf4_table %>%
  dplyr::filter(Group1 == "AncientBritish" & Group2 == "AncientBritish") %>%
  ggplot(aes(x = rasAF, y = abs(rasF4Z))) +
  geom_jitter(width=0.25)
ggsave("analyses/HGDP-exploration/rasF4_1000G_ancientsOnly_tvOnly.pdf")

rasf4_table %>%
  dplyr::filter(Group1 == "AncientBritish" & Group2 == "AncientBritish") %>%
  dplyr::filter(abs(rasF4Z) >= 3) %>%
  dplyr::group_by(rasAF) %>%
  dplyr::summarise(nr_significant = dplyr::n())

# Checking Coverage dependence
ras_table %>%
  dplyr::filter(dataset == "1000G" & rasAF == "All" & Group == "AncientBritish" &
                Right != "YRI2" & Right != "CHB2") %>%
  ggplot(aes(x = Norm, y = RAS, col = Right)) + geom_point() +
  geom_errorbar(aes(ymin = RAS - StdErr, ymax = RAS + StdErr))
ggsave("analyses/HGDP-exploration/coverage_dependence_tvOnly.pdf")
