#!/usr/bin/env Rscript

## Required packages
# library("readr")
# library("rlang")
# library("dplyr", quietly = TRUE)
# library("stringr")
# library("stats")
# library("tibble")
library("tidyr") ## Cannot find a way to use %>% with a package specific call
# library("purrr", quietly = TRUE)
library("foreach")

## Function to get the path of this script. Needed to source Defs.R correctly.
getScriptPath <- function(){
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}

## Source function definitions
source(paste(getScriptPath(),"Defs.R", sep="/"))
# write("Defs loaded.", file=stderr())
## Read command line arguments
args = commandArgs(trailingOnly=TRUE)

## If user asked for helptext, print it and exit
if (args[1] %in% c("-h", "--help", NA)){
  write("Usage: ras_to_rasta.R rasta_setup_file n_ind_per_pop input_rascal_output_file_prefix input_rascal_output_file_suffix
  
  rasta_setup_file                  A file with ind/pops A, B and C to calculate desired rasta(A, B; C, Outgroup)
  
  num_chroms                        The number of chromosomes simulated.
  
  max_ac                            The maximum allele count used in RAS calculations.
  
  n_ind_per_pop                     The number of inidividuals in each simulated population.
  
  input_rascal_output_file_prefix   The filepath to the input files, up to the chromosome number.
  
  input_rascal_output_file_suffix   The suffix of the input files, following the chromosome number.
  
  num_threads                          The number of threads to use for parallelisation
", file=stderr())
quit(status = 1)
}

setup_file <- args[1]
num_chroms <- args[2]
max_ac <- args[3]
n_ind_per_pop <- as.integer(args[4])
prefix <- args[5]
suffix <- args[6]
num_threads <- args[7]


## Make parallel cluster
cl <- parallel::makeCluster(num_threads)
doParallel::registerDoParallel(cl)


## Read in the desired RASTA setups
rasta_setups <- readr::read_tsv(setup_file, col_names = c("A","B","C"), col_types='ccc')

## Load data in one huge tibble
ras_dat <- 1:num_chroms %>% purrr::map_dfr(~load_chrom(prefix, ., suffix))

## Compile a tibble with RASTA/F4 for all tests per AC.
# rasta_dat <- NULL
for (ac in c(2:max_ac, "Outgroup F3")) {
  ## Sum up individual-to-individual sharing to create an individual-to-population sharing. Up to AC 5
  grouped_data <- get_group_sharing(ras_dat, n_ind_per_pop, mode=ac)
  
  ## This will parallelise the operations in the loop and combine the outputs using bind_rows.
  rasta_dat <- foreach::foreach(row=1:nrow(rasta_setups), .combine=dplyr::bind_rows ) %dopar% {
    ## Library loading in parallel chunk needed to use %>%
    library(tidyr)
      calculate_rasta(grouped_data, !!rlang::sym(rasta_setups$A[row]),  !!rlang::sym(rasta_setups$B[row]), !!rlang::sym(rasta_setups$C[row])) %>%
        delete_mj_jackknife(data=., values=rasta, mj=NumSites, weights=NumSites) %>%
        tibble::as_tibble() %>% dplyr::mutate(A=rasta_setups$A[row], B=rasta_setups$B[row], C=rasta_setups$C[row], AlleleCount=ac) %>%
        dplyr::select(A,B,C,theta_J, jack_se, Zscore, AlleleCount)
  }
  ## Save output to file
  save_output("rasta_out.m", ac,".txt", rasta_dat)
}

## Stop the parallel cluster
parallel::stopCluster(cl)
