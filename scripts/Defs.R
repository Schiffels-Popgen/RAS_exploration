## Load chromosome data to memory
load_chrom <- function(prefix,chr,suffix) {
  fn <- paste0(prefix, chr, suffix)
  readr::read_tsv(fn, comment = "#", col_types='ccddddc') %>% dplyr::mutate(chrom = as.character(chr))
}

## Save output into files depending on AC
save_output <- function(prefix, ac, suffix, data) {
  fn <- paste(prefix, stringr::str_replace(ac," ","_"), suffix, sep=".")
  readr::write_tsv(data %>% dplyr::filter(AlleleCount == ac), fn)
}

## Implementation of delete mj jackknifing following Nick Patterson's writeup of it.
## The function works like tidyverse functions, i.e. values, mj and weights should not be quoted/be symbols.
## Returns a list of items :
          # theta_J ## The jackknife mean from all the theta_-x data.
          # jack_se ## The SE of theta_J
          # Zscore ## The Zscore of the theta_J
          # jack.values ## A vector of all the theta_-x values
delete_mj_jackknife <- function(data, values, mj, weights, .print_vals=F) {
  values <- rlang::enquo(values)
  mj <- rlang::enquo(mj)
  weights <- rlang::enquo(weights)
  
  data <- data %>%  dplyr::ungroup() ## Remove grouping from the data to avoid summarise mishaps
  
  ## Calculate the sum of each column as well as theta_hat
  sums <- data %>%  dplyr::summarise(value_sum=sum(!!values),     ## The sum of values
                             mj_sum=sum(!!mj),            ## The sum of mj. used for calculating rates of the values
                             weight_sum=sum(!!weights),   ## The sum of weights across all blocks
                             theta_hat=value_sum/mj_sum   ## Theta hat as a rate of values across all blocks
                             )
  
  ## Add columns with the theta_minus, as well as n_minus_weight 
  data <- data %>% dplyr::mutate(
    value_minus=sums$value_sum-!!values,        ## value total - value of block
    mj_minus=sums$mj_sum-!!mj,                  ## mj total - mj of block
    n_minus_weight=sums$weight_sum-!!weights,   ## weight total - weight of block
    theta_minus=value_minus/mj_minus,           ## The rate calculated from all but the given block.
    hj=sums$weight_sum/!!weights,               ## hj per block for tau calculation
    tau=hj*sums$theta_hat - (hj-1)* theta_minus ## pseudovalue for variance calculation
    )
  ## Use added columns to calculate theta_J
  theta_J <- data %>%
     dplyr::summarise(theta_J=nrow(data) * sums$theta_hat - sum(n_minus_weight*theta_minus/sums$weight_sum)) %>% 
    dplyr::pull(theta_J)
  variance <- data %>%  dplyr::summarise(variance=1/nrow(data) * sum((tau-theta_J)^2/(hj-1))) %>% dplyr::pull(variance)
  jack_se = sqrt(variance)
  if (.print_vals) {
  return(
    list(
          theta_J=theta_J,
          jack_se=jack_se,
          Zscore=theta_J/jack_se,
          jack.values=as.vector(data$theta_minus)
          )
    )
  } else {
    return(
    list(
          theta_J=theta_J,
          jack_se=jack_se,
          Zscore=theta_J/jack_se
          )
    )
  }
}
  
## Function to dplyr::pull sharing between two individuals at a specific AlleleCount.
filter_ras_blocks <- function(ras_dat, A, B, mode = "Total [2,5]") {
  ras_dat %>%
    dplyr::filter(TestPop == A & LeftPop == B & AlleleCount == mode) 
}

## Function to dplyr::pull sharing between two individuals. Sums rare variant allele counts up to the given AC.
filter_ras_blocks_up_to_ac <- function(ras_dat, A, B, mode = 5) {
  ## If mode is Outgroup F3 or total, use mode as is
  if ( any(stringr::str_detect(mode, c("Outgroup", "Total"))) ) {
    ras_dat %>%
    dplyr::filter(TestPop == A & LeftPop == B & AlleleCount == mode) %>%
     dplyr::group_by(TestPop, LeftPop, chrom) %>%
     dplyr::summarise(.groups="keep", RAS=sum(RAS), NumSites=dplyr::first(NumSites), `RAS/site`=sum(RAS)/NumSites, Error=dplyr::first(Error), AlleleCount=mode) %>%
    dplyr::arrange(chrom)
  ## Else, sum up to given AC
  } else {
  ras_dat %>%
    dplyr::filter(TestPop == A & LeftPop == B & AlleleCount %in% c(2:mode)) %>%
     dplyr::group_by(TestPop, LeftPop, chrom) %>%
     dplyr::summarise(.groups="keep", RAS=sum(RAS), NumSites=dplyr::first(NumSites), `RAS/site`=sum(RAS)/NumSites, Error=dplyr::first(Error), AlleleCount=mode) %>%
    dplyr::arrange(chrom)
  }
}

## A function to fetch the within-population/individual sharing value up to a given AC
fetch_intra_RAS_value <- function(ras_dat, ind_chr_df, mode = 5) {
  x <- NULL
  for (row in 1:nrow(ind_chr_df)) {
    ind=ind_chr_df[[row,1]]
    chr=ind_chr_df[[row,2]]
    # print(ind)
    # print(chr)
    x <- x %>% dplyr::bind_rows(filter_ras_blocks_up_to_ac(ras_dat, A=ind, B=ind, mode=mode) %>%
    dplyr::filter(chrom == chr)) 
  }
  return(x %>% dplyr::pull(RAS))
}

## Function to merge individual sharings into an uncorrected group sharing
merge_individual_sharing <- function(ras_dat, n_ind_per_pop = 10, mode = 5) {
  groups <- c(paste0("Pop",c(0:8)))
  x <- c(paste0("ind",c(0:(n_ind_per_pop*9-1))))
  y <- rep(groups, each=n_ind_per_pop)
  groupings <- tibble::tibble(ind=x,pop=y)
  ## Allow querying for either Outgroup F3/Total RAS from RAS table
  if ( any(stringr::str_detect(mode, c("Outgroup", "Total"))) ) {
    ac_query <- mode
  ## Or sum ras up to a given ac value
  } else {
    ac_query <- c(2:mode)
  }
  ras_dat %>%
    dplyr::filter(AlleleCount %in% ac_query) %>%
     dplyr::group_by(TestPop, LeftPop, chrom) %>%
    ## Sum the RAS up to the given AF
     dplyr::summarise(.groups="keep", RAS=sum(RAS), NumSites=dplyr::first(NumSites), `RAS/site`=sum(RAS)/NumSites, 
              AlleleCount=mode) %>%
    ## Add popgroup column
    dplyr::left_join(groupings, by=c("TestPop"="ind")) %>%
     dplyr::group_by(LeftPop, pop, chrom) %>% 
    ## Get sharing of each individual with each pop
     dplyr::summarise(.groups="keep", RAS=sum(RAS), NumSites=dplyr::first(NumSites), AlleleCount=dplyr::first(AlleleCount)) %>%
    dplyr::left_join(groupings, by=c("LeftPop"="ind")) %>%
    dplyr::rename(TestPop=pop.x, truePopLeft=pop.y)
}

## Takes the merged individual sharing and corrects it to exclude the within individual sharing of the LeftPop individual before noramlising by n of the population.
get_group_sharing <- function(ras_dat, n_ind_per_pop = 10, mode = 5) {
  ## If the individual belongs to the population it is checked against, subtract within individual sharing from total.
  f1<- merge_individual_sharing(ras_dat, n_ind_per_pop, mode) %>% 
    dplyr::filter(truePopLeft == TestPop) %>%
     dplyr::ungroup() %>%
    ## RAS needs to be divided by the number of individuals in the pop (excluding the test indA)
    dplyr::mutate(RAS=(RAS-fetch_intra_RAS_value(ras_dat, tibble::tibble(x=LeftPop, y=chrom), mode= mode))/(n_ind_per_pop-1))
  ## Keep sharing with other populations unchanged
  f2 <- merge_individual_sharing(ras_dat, n_ind_per_pop, mode) %>% 
    dplyr::filter(truePopLeft != TestPop) %>%
    ## RAS needs to be divided by the number of individuals in the pop
    dplyr::mutate(RAS=RAS/n_ind_per_pop) %>%
     dplyr::ungroup()
  dat <- dplyr::bind_rows(f1,f2) %>%  dplyr::group_by(LeftPop, TestPop, chrom)
  dat
}

## Playground function to sum sharing across chromosomes and print the sd of the distribution.
get_total <- function(ras_dat, A, B, mode = "Total [2,5]") {
  filter_ras_blocks(ras_dat, A, B, mode) %>%
     dplyr::summarise(TestPop=A, LeftPop=B, RAS=sum(RAS), NumSites=sum(NumSites), Error=stats::sd(`RAS/site`)/sqrt(20), `RAS/site`=RAS/NumSites, AlleleCount=mode, chrom="all")
}

## Function that takes the output of get_group_sharing and calculates RASTA(A,C) - RAS(B,C).
## Works like tidyverse functions, i.e. A, B and C should not be quoted/be symbols.
## 'A' should be an individual, while B and C populations.
calculate_rasta <- function(data, A, B, C) {
  # print(A)
  A <- rlang::enquo(A)
  B <- rlang::enquo(B)
  C <- rlang::enquo(C)

  ## Check that A is an individual and B and C are populations, else throw an error
  stopifnot(stringr::str_detect(rlang::quo_name(A),"ind"))
  stopifnot(stringr::str_detect(c(rlang::quo_name(B),rlang::quo_name(C)),"Pop"))
  
  ## First, take only the RAS between indA and popB
  ras_A_C <- data %>%
    dplyr::filter(LeftPop == rlang::quo_name(A), TestPop == rlang::quo_name(C))
  ## To calculate RAS between popB and C, sum RAS of all individuals in PopB that are NOT indA
  ras_B_C <- data %>%
    dplyr::filter(LeftPop != rlang::quo_name(A), TestPop == rlang::quo_name(C), truePopLeft == rlang::quo_name(B)) %>%
     dplyr::group_by(truePopLeft, TestPop, chrom) %>%
     dplyr::summarise(
      .groups="keep",
      LeftPop=rlang::quo_name(B), 
      TestPop=dplyr::first(TestPop), 
      chrom=dplyr::first(chrom), 
      ## RAS needs to be divided by the number of individuals being summed as it is the AVERAGE sharing an individual from the population.
      RAS=sum(RAS)/dplyr::n(), 
      NumSites=dplyr::first(NumSites), 
      AlleleCount=dplyr::first(AlleleCount)
      )
  ## Merge the two tibbles, and calculate the difference of RAS
  rasta_dat <- dplyr::bind_rows(ras_A_C,ras_B_C) %>%
    tidyr::pivot_wider(id_cols=-truePopLeft, names_from = LeftPop, values_from= RAS) %>% 
    dplyr::mutate(rasta=((!!A-!!B)/`NumSites`))
  rasta_dat
}