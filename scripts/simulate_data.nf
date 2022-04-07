#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  =========================================
  rare_variant_data_simulation_analysis v${workflow.manifest.version}
  =========================================
  Usage:

  The typical command for running the pipeline on sdag is as follows:

  nextflow run pipeline.nf -profile cobra,mpcdf,conda,var_sim --four_mN 1 --chrom_length 1e6 --n_ind_per_pop 20 --max_ras_ac 5 --knn 5

  Mandatory arguments:
      -profile [str]          Institution or personal hardware config to use (e.g. standard, docker, singularity, conda, aws). Ask your system admin if unsure, or check documentation.

      --four_mN [float]           The scaled migration rate between non-diagonal neighbour populations used in the simulation.

      --chrom_length [float]  The length of the simulated chromosomes.

      --n_ind_per_pop [int]   The sample size of each of the 9 populations. All populations have an equal population and sample size.
  """.stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
/* --                SET UP CONFIGURATION VARIABLES                       -- */
///////////////////////////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Small console separator to make it easier to read errors after launch
println ""

//////////////////////////////////////////////////////////////
/* --          Create chromosome channel                 -- */
//////////////////////////////////////////////////////////////
/*chromosome_names = Channel.of(1..20)*/

process msprime{
  
  tag "m${params.four_mN}_chr${chrom_name}_l${params.chrom_length}"
  publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${params.chrom_length}/${params.four_mN}", mode: 'copy'
  memory '8GB'
  
  input:
  val chrom_name from Channel.of(1..20)
  
  output:
  tuple val(chrom_name), val("all"), path("all_vars_m${params.four_mN}_chr${chrom_name}.geno"), path("all_vars_m${params.four_mN}_chr${chrom_name}.snp"), path("all_vars_m${params.four_mN}_chr${chrom_name}.ind"), path("all_vars_m${params.four_mN}_chr${chrom_name}.indEach") into (ch_all_vars_datasets, ch_all_vars_for_trident)
  tuple val(chrom_name), val("common"), path("common_vars_m${params.four_mN}_chr${chrom_name}.geno"), path("common_vars_m${params.four_mN}_chr${chrom_name}.snp"), path("common_vars_m${params.four_mN}_chr${chrom_name}.ind"), path("common_vars_m${params.four_mN}_chr${chrom_name}.indEach") into (ch_common_vars_datasets, ch_for_1240k_input_geno, ch_for_1240k_input_snp, ch_for_1240k_input_ind, ch_for_1240k_input_indEach)
  tuple val(chrom_name), val("rare"), path("all_vars_m${params.four_mN}_chr${chrom_name}.freqsum.gz") into (ch_freqsum_dataset, ch_for_linecount)
/*  tuple val(chrom_name), path("all_vars_m${params.four_mN}_chr${chrom_name}.geno") into ch_all_vars_geno_for_f3
  tuple val(chrom_name), path("all_vars_m${params.four_mN}_chr${chrom_name}.snp") into ch_all_vars_snp_for_f3
  tuple val(chrom_name), path("all_vars_m${params.four_mN}_chr${chrom_name}.ind") into ch_all_vars_ind_for_f3
  tuple val(chrom_name), path("all_vars_m${params.four_mN}_chr${chrom_name}.indEach") into ch_all_vars_indEach_for_f3
  tuple val(chrom_name), path("common_vars_m${params.four_mN}_chr${chrom_name}.geno") into (ch_common_vars_geno_for_f3, ch_common_vars_geno_for_1240k )
  tuple val(chrom_name), path("common_vars_m${params.four_mN}_chr${chrom_name}.snp") into (ch_common_vars_snp_for_f3, ch_common_vars_snp_for_1240k )
  tuple val(chrom_name), path("common_vars_m${params.four_mN}_chr${chrom_name}.ind") into (ch_common_vars_ind_for_f3, ch_common_vars_ind_for_1240k )
  tuple val(chrom_name), path("common_vars_m${params.four_mN}_chr${chrom_name}.indEach") into (ch_common_vars_indEach_for_f3, ch_common_vars_indEach_for_1240k )*/
  path("variant_counts.m${params.four_mN}_chr${chrom_name}.txt") 
  
  script:
  """
  ${baseDir}/simulation.py -c ${params.chrom_length} -n ${chrom_name} -s ${params.n_ind_per_pop} -m ${params.four_mN} #-o /projects1/MICROSCOPE/rarevar_sim_study/data/
  ${baseDir}/add_ref_and_sort_eigenstrat.sh ${chrom_name} ${params.four_mN} .
  """
}

/*
Collect and filter the 1240k making input into 4 subchannels
*/
ch_common_vars_geno_for_1240k=ch_for_1240k_input_geno.map{ it[2] }
ch_common_vars_snp_for_1240k=ch_for_1240k_input_snp.map{ it[3] }
ch_common_vars_ind_for_1240k=ch_for_1240k_input_ind.map{ it[4] }
ch_common_vars_indEach_for_1240k=ch_for_1240k_input_indEach.map{ it[5] }


process make_1240k{

  tag "m${params.four_mN}_chr${chrom_name}_l${params.chrom_length}"
  publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${params.chrom_length}/${params.four_mN}", mode: 'copy'
  memory '1GB'

  input:
  path genos from ch_common_vars_geno_for_1240k.collect()
  path snps from ch_common_vars_snp_for_1240k.collect()
  path inds from ch_common_vars_ind_for_1240k.collect()
  path indEachs from ch_common_vars_indEach_for_1240k.collect()

  output:
  tuple val("1:20"),val("1240k"),path("twelve_forty_vars_m${params.four_mN}.geno"),path("twelve_forty_vars_m${params.four_mN}.snp"),path("twelve_forty_vars_m${params.four_mN}.ind"),path("twelve_forty_vars_m${params.four_mN}.indEach") into ch_twelve_forty_vars_datasets

  script:
  """
  ## First concatenate the snps and genos by chromosome
  for chrom_name in {1..20}; do
    cat common_vars_m${params.four_mN}_chr\${chrom_name}.geno  >> concat_genos
    cat common_vars_m${params.four_mN}_chr\${chrom_name}.snp  >> concat_snps
  done
  
  ## Then paste side by side, use shuf to randomly subsample 1200K sites, sort by original line numbers, and remove leading spaces (added by cat -n)
  paste -d " " concat_snps concat_genos | cat -n | shuf -n 1200000 | sort -nk1 | sed -e 's/^[ ]*//' >banana.txt
  
  ## Extract the .snp and .geno file from the subsample. 
  cut -f 1 -d " " banana.txt | cut -f 2- >twelve_forty_vars_m${params.four_mN}.snp
  cut -f 2 -d " " banana.txt >twelve_forty_vars_m${params.four_mN}.geno
  
  ## Copy the original ind files 
  cp ${inds[0]} twelve_forty_vars_m${params.four_mN}.ind
  cp ${indEachs[0]} twelve_forty_vars_m${params.four_mN}.indEach
  """
}

