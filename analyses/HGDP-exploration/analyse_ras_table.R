library(magrittr)

ras_table <- readr::read_delim("analyses/HGDP-exploration/AncientBritish_HGDP_ras.table.txt", 
                               delim = "\t", escape_double = FALSE, trim_ws = TRUE)

ras_table %>% dplyr::filter(Right=="Basque2")

