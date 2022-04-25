library(ggplot2)
source("Script/Functions/load_data.R")

# ras_table <- read_ras_tables()
# readr::write_rds(ras_table, "Data/combined_data.rds")
ras_table <- readr::read_rds("Data/combined_data.rds")

ras_table %>% dplyr::filter(Right == "CEU2" & dataset == "1000G" &
                             TF == FALSE & rasAF == "01" & tvOnly == TRUE &
                             mapMasked == TRUE)


ras_table %>% dplyr::filter(
  dataset == "1000G",
  TF == FALSE,
  mapMasked == TRUE,
  tvOnly == TRUE,
  Right == "FIN2",
  Left %in% c("<12881A>", "<12884A>")
) %>% dplyr::select(Left, rasAF, RAS, StdErr) %>%
  tidyr::pivot_wider(names_from = Left, values_from = c(RAS, StdErr)) %>%
  mutate(
    rasF3 = `RAS_<12881A>` - `RAS_<12884A>`,
    
  )

# Snippets

# Checking dependency on Coverage
ras_table %>%
  dplyr::filter(dataset == "1000G" & rasAF == "All" & Group == "AncientBritish" &
                  Right != "YRI2" & Right != "CHB2") %>%
  ggplot(aes(x = Norm, y = RAS, col = Right)) + geom_point() +
  geom_errorbar(aes(ymin = RAS - StdErr, ymax = RAS + StdErr))
ggsave("analyses/HGDP-exploration/coverage_dependence.pdf")

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
  dplyr::filter(!(Group1 %in% c("AncientBritish", "YRI", "CHB", "PEL")),
                !(Group2 %in% c("AncientBritish", "YRI", "CHB", "PEL"))) %>%
  ggplot(aes(x = rasAF, y = abs(rasF4Z))) +
  geom_jitter(width=0.25)
ggsave("analyses/HGDP-exploration/rasF4_1000G_modernOnly.pdf")

rasf4_table %>%
  dplyr::filter(!(Group1 %in% c("AncientBritish", "YRI", "CHB", "PEL")),
                !(Group2 %in% c("AncientBritish", "YRI", "CHB", "PEL")),
                abs(rasF4Z) > 3) %>%
  dplyr::group_by(rasAF) %>% dplyr::summarise(dplyr::n())

rasf4_table %>%
  dplyr::filter(Group1 == "AncientBritish" & Group2 == "AncientBritish") %>%
  ggplot(aes(x = rasAF, y = abs(rasF4Z))) +
  geom_jitter(width=0.25)
ggsave("analyses/HGDP-exploration/rasF4_1000G_ancientsOnly.pdf")

rasf4_table %>%
  dplyr::filter(Group1 == "AncientBritish" & Group2 == "AncientBritish") %>%
  dplyr::filter(abs(rasF4Z) >= 3) %>%
  dplyr::group_by(rasAF) %>%
  dplyr::summarise(nr_significant = dplyr::n())
