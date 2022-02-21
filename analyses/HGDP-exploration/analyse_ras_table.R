library(magrittr)
library(ggplot2)

ras_table <- readr::read_tsv("analyses/HGDP-exploration/AncientBritish_HGDP_ras.table.txt")

ras_table %>% dplyr::filter(Cumulative==TRUE, k==5) %>% ggplot() +
  geom_errorbar(aes(x=Left, y=RAS, ymin=RAS-StdErr, ymax=RAS+StdErr, col=Right)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

