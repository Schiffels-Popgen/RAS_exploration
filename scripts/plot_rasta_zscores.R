#!/usr/bin/env Rscript

library("tidyr") ## Cannot find a way to use %>% with a package specific call

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

## Read command line arguments
args = commandArgs(trailingOnly=TRUE)

## If user asked for helptext, print it and exit
if (args[1] %in% c("-h", "--help", NA)){
  write("Usage: ras_to_rasta.R rasta_setup_file n_ind_per_pop input_rascal_output_file_prefix input_rascal_output_file_suffix
  
  chrom_length                      The chromosome length used in the simulations.
  
  max_ac                            The maximum allele count used in RAS calculations.
  
  n_ind_per_pop                     The number of inidividuals in each simulated population.
  
  input_rascal_output_file_prefix   The filepath to the input files, up to but not including the allele count.
  
  input_rascal_output_file_suffix   The suffix of the input files, following the allele count.
", file=stderr())
quit(status = 1)
}

chrom_length <- args[1]
max_ac <- args[2]
n_ind_per_pop <- as.integer(args[3])
prefix <- args[4]
suffix <- args[5]

rasta_data <- c(2:max_ac,"Outgroup_F3") %>% purrr::map_dfr(~load_rasta(prefix, ., suffix))

rasta_data <- categorise_rasta(rasta_data, n_ind_per_pop)

ggplot2::ggplot(rasta_data, ggplot2::aes(x=Type, y=abs(Zscore))) +
  ggplot2::scale_x_discrete() + ## Needed to place geom_hline below geom_point without getting errors
  ggplot2::geom_vline(xintercept=seq(0,8,1)+.5,color="gray40",alpha=0.5, size=1) +  ## Vertical line separators between the different rasta types.
  ggplot2::geom_hline(yintercept=3,color="darkred", alpha=0.4, size=1) +  ## Marks the |Z| = 3 line
  ggplot2::geom_point(ggplot2::aes(x=Type, y=abs(Zscore)), position=ggplot2::position_jitter(w=0.4,h=0), alpha=0.6, shape=21, col="black", fill="grey40", stroke=0.7, size=2) +
  ggplot2::theme_bw() +
  ggplot2::xlab('') + 
  ggplot2::ylab('| Z-Score |') +
  ggplot2::theme(panel.grid.minor.x=ggplot2::element_line(color="grey50"), axis.text.x=ggplot2::element_text(angle=90, hjust = 1), axis.ticks.x = ggplot2::element_blank()) +
  ggplot2::facet_wrap(~AlleleCount, scales='free_x', ncol=3)

ggplot2::ggsave(paste0("Rasta_ZScore_plot_l",chrom_length,".pdf"), unit="cm", width=30, height=20)
