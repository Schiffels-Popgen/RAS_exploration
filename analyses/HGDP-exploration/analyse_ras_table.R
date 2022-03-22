library(magrittr)
library(ggplot2)

ras_table <- readr::read_tsv("analyses/HGDP-exploration/AncientBritish_HGDP_ras.table.txt")

# Complete view
ras_table %>% dplyr::filter(Cumulative==TRUE, k==5) %>% ggplot() +
  geom_errorbar(aes(x=Left, y=RAS, ymin=RAS-StdErr, ymax=RAS+StdErr, col=Right)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# One left vs multiple right
ras_table %>% dplyr::filter(Cumulative & k==2 & Left=="<12881A>") %>% ggplot() +
  geom_errorbar(aes(x=Right, y=RAS, ymin=RAS-StdErr, ymax=RAS+StdErr)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# one right vs multiple left
ras_table %>% dplyr::filter(Cumulative & k==2 & Right=="French2") %>% ggplot() +
  geom_errorbar(aes(x=Left, y=RAS, ymin=RAS-StdErr, ymax=RAS+StdErr)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# vs allele freq
ras_table %>% dplyr::filter(grepl("1????A", Left) & !Cumulative & Right=="French2") %>% ggplot() +
  geom_errorbar(aes(x=k, y=RAS, ymin=RAS-StdErr, ymax=RAS+StdErr, col=Left)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ras_table %>% dplyr::filter(!Cumulative & Left=="<12881A>") %>% ggplot() +
  geom_errorbar(aes(x=k, y=RAS, ymin=RAS-StdErr, ymax=RAS+StdErr, col=Right)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
