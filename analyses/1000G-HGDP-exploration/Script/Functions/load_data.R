library(magrittr)

file_pattern <- "AncientBritish_(HGDP|1000G)_(1240K_){0,1}ras([0-9A-Za-z]+)(_TVonly){0,1}(_mapMasked){0,1}.table.txt"

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
  "<HG02425>", "PEL",
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
  m <- stringr::str_match(fn, file_pattern)
  dataset <- m[2]
  TF <- if(is.na(m[3])) FALSE else TRUE
  af <- m[4]
  tvOnly <- if(is.na(m[5])) FALSE else TRUE
  mapMasked <- if(is.na(m[6])) FALSE else TRUE
  readr::read_tsv(full_filename, col_types="ccidd") %>%
    dplyr::mutate(dataset = dataset, TF = TF, rasAF = af, tvOnly = tvOnly, mapMasked = mapMasked) %>%
    dplyr::left_join(popNames)
}

read_ras_tables <- function(dir = "Data") {
  list.files(dir, pattern = file_pattern, full.names=TRUE) %>%
    purrr::map_dfr(~read_ras_table(.))
}
