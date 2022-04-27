library(ggplot2)
source("Script/Functions/load_data.R")

# ras_table <- read_ras_tables()
# readr::write_rds(ras_table, "Data/combined_data.rds")
ras_table <- readr::read_rds("Data/combined_data.rds")

# 2D plot
ras_table %>%
  dplyr::filter(dataset == "1000G", TF == FALSE, rasAF == "Common", tvOnly == TRUE, mapMasked == TRUE) %>%
  dplyr::select(Left, Right, RAS, StdErr, Group) %>%
  dplyr::filter(!(Group %in% c("AncientBritish"))) %>%
  tidyr::pivot_wider(names_from = Right, values_from = c(RAS, StdErr)) %>%
  ggplot(aes(x = RAS_IBS2, y = RAS_FIN2, color = Group)) +
    geom_errorbar(aes(ymin = RAS_FIN2 - StdErr_FIN2, ymax = RAS_FIN2 + StdErr_FIN2)) +
    geom_errorbarh(aes(xmin = RAS_IBS2 - StdErr_IBS2, xmax = RAS_IBS2 + StdErr_IBS2))

# 1D
ras_table %>%
  dplyr::filter(dataset == "1000G", TF == FALSE, rasAF == "01", tvOnly == TRUE, mapMasked == TRUE) %>%
  dplyr::select(Left, Right, RAS, StdErr, Group) %>%
  # dplyr::filter(!(Group %in% c("AncientBritish", "YRI"))) %>%
  dplyr::filter(Group == "AncientBritish") %>%
  tidyr::pivot_wider(names_from = Right, values_from = c(RAS, StdErr)) %>%
  ggplot(aes(x = Left, y = RAS_FIN2 - RAS_IBS2, color = Group)) +
  geom_errorbar(aes(ymin = RAS_FIN2 - StdErr_FIN2, ymax = RAS_FIN2 + StdErr_FIN2))


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
ras_table_short <- ras_table %>% dplyr::filter(
  dataset == "1000G",
  TF == FALSE,
  tvOnly == FALSE,
  mapMasked == FALSE
) %>% dplyr::select(Left, Right, rasAF, RAS, StdErr, Group)

inds <- ras_table_short %>% dplyr::pull(Left) %>% unique()
rasAF <- ras_table %>% dplyr::pull(rasAF) %>% unique()
rasf4_table <- tidyr::expand_grid(Left1 = inds, Left2 = inds, rasAF = rasAF) %>%
  dplyr::filter(Left1 != Left2) %>%
  dplyr::left_join(ras_table_short, by = c("Left1" = "Left", "rasAF" = "rasAF")) %>%
  dplyr::rename(RAS1 = RAS, StdErr1 = StdErr, Group1 = Group) %>%
  dplyr::left_join(ras_table_short, by = c("Left2" = "Left", "rasAF" = "rasAF", "Right" = "Right")) %>%
  dplyr::rename(RAS2 = RAS, StdErr2 = StdErr, Group2 = Group) %>%
  dplyr::mutate(
    f4name = paste0("RAS-F4(", Left1, ", ", Left2, ", ", Right, ", Outgroup)"),
    rasF4 = RAS1 - RAS2,
    rasF4err = sqrt(StdErr1^2 + StdErr2^2),
    rasF4Z = rasF4 / rasF4err
  )

excl <- c("AncientBritish", "YRI", "CHB", "PEL")
rasf4_table %>%
  dplyr::filter(!(Group1 %in% excl),
                !(Group2 %in% excl),
                Group1 != substr(Right, 1, 3),
                Group2 != substr(Right, 1, 3)) %>%
  ggplot(aes(x = rasAF, y = abs(rasF4Z), col = (Group1 == Group2))) +
  geom_jitter(width=0.25)
ggsave("Output/rasF4_1000G_modernOnly_groupCol.pdf")

incl <- c("<12880A>", "<12881A>", "<12883A>", "<12884A>", "<12885A>", 
          "<15558A>", "<15569A>", "<15570A>", "<15577A>")
ia <- c("<12880A>", "<12884A>", "<15579A>", "<M1489>")
rasf4_table %>%
  dplyr::filter(Left1 %in% incl,
                Left2 %in% incl) %>%
  ggplot(aes(x = rasAF, y = abs(rasF4Z), col = Left1)) +
  geom_jitter(width=0.25)
ggsave("Output/rasF4_1000G_ancientsOnly.pdf")

ras_table_short %>%
  dplyr::filter(Right == "CEU2",
                rasAF == "01",
                !(Group %in% c("YRI", "PEL", "CHB")),
                Left %in% c("1288") %>%
  ggplot(aes(x = Left, y = RAS, col = Group)) +
    geom_errorbar(aes(ymin = RAS - StdErr, ymax = RAS + StdErr))


