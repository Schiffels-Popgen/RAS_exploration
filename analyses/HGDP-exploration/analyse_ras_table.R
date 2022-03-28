library(magrittr)
library(ggplot2)

pat <- "AncientBritish_(HGDP|1000G)_ras([0-9A-Za-z]+).table.txt"

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

#plot1 <-
ras_table %>% dplyr::filter(dataset == "1000G" & rasAF == "01" & !(Right %in% c("CHB2", "YRI2")) &
                              !(Group %in% c("YRI", "CHB", "PEL"))) %>%
  dplyr::select(Left, Group, Right, RAS, StdErr) %>%
  ggplot(aes(x = Left, y = RAS, col = Group)) + geom_point() +
  geom_errorbar(aes(ymin = RAS - StdErr, ymax = RAS + StdErr)) +
  facet_wrap(~Right) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("sanalyses/HGDP-exploration/scatter_plot_1000G_ras01.pdf", plot1)

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
  fn <- paste0("analyses/HGDP-exploration/scatter_plot_1000G_ras", freq, "_noFin.pdf")
  ggsave(fn, plot)
}


# Playing with PCA and Umap, but isn't really useful
leftRight <- ras_table %>% dplyr::filter(dataset == "1000G" & rasAF == "01" & !(Right %in% c("CHB2", "YRI2")) &
                              !(Group %in% c("YRI", "CHB", "PEL", "FIN"))) %>%
  dplyr::select(Left, Group, Right, RAS) %>%
  tidyr::pivot_wider(names_from = Right, values_from = RAS)
ind_names <- leftRight$Left
group_names <- leftRight$Group
mat <- dplyr::select(leftRight, -c(Left, Group)) %>% as.matrix() %>% t()
pca <- prcomp(mat)
pca_dat <- tibble::as.tibble(pca$rotation) %>%
  dplyr::mutate(Ind = ind_names, Group = group_names)

ggplot(pca_dat, aes(x = PC1, y = PC2, col = Group)) + geom_point()

umap_obj <- umap::umap(t(mat))
umap_dat <- tibble::as.tibble(umap_obj$layout) %>%
  dplyr::mutate(Ind = ind_names, Group = group_names)
  
ggplot(umap_dat, aes(x = V1, y = V2, col = Group)) + geom_point()


# Checking F4-values


