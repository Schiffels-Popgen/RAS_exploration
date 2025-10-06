#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  =========================================
  rare_variant_data_simulation_analysis v${workflow.manifest.version}
  =========================================
  Usage:

  Runs a BACKWARDS simulation in msprime with a change in migration rate at time t_m_change, for 9 populations with a specified demographic history. 
  Then adds a reference individual to the data, sorts the data into Eigenstrat format, and creates four datasets: all variants, rare variants with allele count up to 10, common variants with MAF>0.05, and a random subset of 1.2 million variants (to mimic the 1240K SNP set).
  Finally, creates Poseidon packages for each of the three datasets.
  The typical command for running the pipeline on sdag is as follows:

  nextflow run simulate_data_chained.nf -profile eva,grace,conda,var_sim --four_mN 1.0 --four_mN2 10.0 --chrom_length 1e6 --n_ind_per_pop 20

  Mandatory arguments:
      -profile [str]          Institution or personal hardware config to use (e.g. standard, docker, singularity, conda, aws). Ask your system admin if unsure, or check documentation.

      --four_mN [float]       The scaled per-generation migration rate for the proximal part of the simulation (t=0 to t=t_m_change).

      --four_mN2 [float]      The scaled per-generation migration rate for the distal part of the simulation (t=t_m_change to end).

      --chrom_length [float]  The length of the simulated chromosomes.

      --t_m_change [int]      The time (in generations) at which the migration rate changes from four_mN to four_mN_2.

      --n_ind_per_pop [int]   The sample size of each of the 9 populations. All populations have an equal population and sample size.
  """.stripIndent()
}

nextflow.enable.dsl=1 // Force DSL1 syntax

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
  
  tag "n${params.n_ind_per_pop}_t${params.t_m_change}_m${params.four_mN}_M${params.four_mN2}_chr${chrom_name}_l${params.chrom_length}"
  // publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${params.chrom_length}/t${params.t_m_change}/${params.four_mN}_${params.four_mN2}", mode: 'copy'
  memory '8GB'
  // executor 'local'

  input:
  val chrom_name from Channel.of(1..20)
  
  output:
  path("all_varsm${params.four_mN}_M${params.four_mN2}_chr${chrom_name}.geno") into ch_all_vars_geno
  path("all_varsm${params.four_mN}_M${params.four_mN2}_chr${chrom_name}.snp") into ch_all_vars_snp
  path("all_varsm${params.four_mN}_M${params.four_mN2}_chr${chrom_name}.ind") into ch_all_vars_ind
  path("common_varsm${params.four_mN}_M${params.four_mN2}_chr${chrom_name}.geno") into (ch_common_vars_geno, ch_common_vars_geno_for_1240k )
  path("common_varsm${params.four_mN}_M${params.four_mN2}_chr${chrom_name}.snp") into (ch_common_vars_snp, ch_common_vars_snp_for_1240k )
  path("common_varsm${params.four_mN}_M${params.four_mN2}_chr${chrom_name}.ind") into (ch_common_vars_ind, ch_common_vars_ind_for_1240k )
  path("variant_counts.m${params.four_mN}_M${params.four_mN2}_chr${chrom_name}.txt") 
  
  script:
  """
  ${baseDir}/chained_simulation.py -c ${params.chrom_length} -n ${chrom_name} -s ${params.n_ind_per_pop} -t ${params.t_m_change} -m ${params.four_mN} -m2 ${params.four_mN2}
  ${baseDir}/add_ref_and_sort_eigenstrat.sh ${chrom_name} ${params.four_mN}_M${params.four_mN2} .
  """
}

process make_1240k{
  tag "n${params.n_ind_per_pop}_t${params.t_m_change}_m${params.four_mN}_M${params.four_mN2}_l${params.chrom_length}"
  publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${params.chrom_length}/t${params.t_m_change}/${params.four_mN}_${params.four_mN2}", mode: 'copy'
  memory '1GB'
  executor 'local'

  input:
  path genos from ch_common_vars_geno_for_1240k.collect()
  path snps from ch_common_vars_snp_for_1240k.collect()
  path inds from ch_common_vars_ind_for_1240k.collect()

  output:
  tuple val("1:20"),val("twelve_forty_vars"),path("twelve_forty_vars.geno"),path("twelve_forty_vars.snp"),path("twelve_forty_vars.ind") into ch_twelve_forty_vars_datasets

  script:
  """
  ## First concatenate the snps and genos by chromosome
  for chrom_name in {1..20}; do
    cat common_varsm${params.four_mN}_M${params.four_mN2}_chr\${chrom_name}.geno  >> concat_genos
    cat common_varsm${params.four_mN}_M${params.four_mN2}_chr\${chrom_name}.snp  >> concat_snps
  done
  
  ## Then paste side by side, use shuf to randomly subsample 1200K sites, sort by original line numbers, and remove leading spaces (added by cat -n)
  paste -d " " concat_snps concat_genos | cat -n | shuf -n 1200000 | sort -nk1 | sed -e 's/^[ ]*//' >banana.txt
  
  ## Extract the .snp and .geno file from the subsample. 
  cut -f 1 -d " " banana.txt | cut -f 2- >twelve_forty_vars.snp
  cut -f 2 -d " " banana.txt >twelve_forty_vars.geno
  
  ## Copy the original ind files 
  cp ${inds[0]} twelve_forty_vars.ind
  """
}

process merge_chroms_common_vars {
  tag "n${params.n_ind_per_pop}_t${params.t_m_change}_m${params.four_mN}_M${params.four_mN2}_l${params.chrom_length}"
  publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${params.chrom_length}/t${params.t_m_change}/${params.four_mN}_${params.four_mN2}/", mode: 'copy'
  memory '50MB'
  cpus 1
  executor 'local'

  input:
  path genos from ch_common_vars_geno.collect()
  path snps from ch_common_vars_snp.collect()
  path inds from ch_common_vars_ind.collect()

  output:
  tuple val("1:20"),val("common_vars"),path("common_vars.geno"),path("common_vars.snp"),path("common_vars.ind") into ch_common_vars_datasets

  script:
  """
  ## First concatenate the snps and genos by chromosome
  for chrom_name in {1..20}; do
    cat common_varsm${params.four_mN}_M${params.four_mN2}_chr\${chrom_name}.geno  >> common_vars.geno
    cat common_varsm${params.four_mN}_M${params.four_mN2}_chr\${chrom_name}.snp  >> common_vars.snp
  done

  ## Copy the original ind files 
  cp ${inds[0]} common_vars.ind
  """
}


process merge_chroms_all_vars {
  tag "n${params.n_ind_per_pop}_t${params.t_m_change}_m${params.four_mN}_M${params.four_mN2}_l${params.chrom_length}"
  publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${params.chrom_length}/t${params.t_m_change}/${params.four_mN}_${params.four_mN2}/", mode: 'copy'
  memory '50MB'
  cpus 1
  executor 'local'

  input:
  path genos from ch_all_vars_geno.collect()
  path snps from ch_all_vars_snp.collect()
  path inds from ch_all_vars_ind.collect()

  output:
  tuple val("1:20"),val("all_vars"),path("all_vars.geno"),path("all_vars.snp"),path("all_vars.ind") into ch_all_vars_datasets

  script:
  """
  ## First concatenate the snps and genos by chromosome
  for chrom_name in {1..20}; do
    cat all_varsm${params.four_mN}_M${params.four_mN2}_chr\${chrom_name}.geno  >> all_vars.geno
    cat all_varsm${params.four_mN}_M${params.four_mN2}_chr\${chrom_name}.snp  >> all_vars.snp
  done

  ## Copy the original ind files 
  cp ${inds[0]} all_vars.ind
  """
}

ch_datasets=ch_all_vars_datasets
    .mix(ch_common_vars_datasets, ch_twelve_forty_vars_datasets)
    .dump(tag:"datasets")

process create_poseidon_packages {
  tag "${variant_set}_n${params.n_ind_per_pop}_t${params.t_m_change}_m${params.four_mN}_M${params.four_mN2}_l${params.chrom_length}"
  publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${params.chrom_length}/t${params.t_m_change}/${params.four_mN}_${params.four_mN2}/poseidon", mode: 'copy'
  memory '50MB'
  cpus 1
  executor 'local'

  input:
  tuple chroms, variant_set, path(geno), path(snp), path(ind) from ch_datasets

  output:
  tuple path("${variant_set}/${variant_set}.bed"), path("${variant_set}/${variant_set}.bim"), path("${variant_set}/${variant_set}.fam"), path("${variant_set}/${variant_set}.geno"), path("${variant_set}/${variant_set}.snp"), path("${variant_set}/${variant_set}.ind"), path("${variant_set}/${variant_set}.janno"), path("${variant_set}/${variant_set}.bib"), path("${variant_set}/POSEIDON.yml")

  script:
  """
  ## trident creates the package within the work directory, and nextflow is responsible for putting it in the data dir.
  ${params.poseidon_exec_dir}/trident init --inFormat EIGENSTRAT \
      --snpSet Other \
      --genoFile ${variant_set}.geno \
      --snpFile ${variant_set}.snp \
      --indFile ${variant_set}.ind \
      -o ${variant_set} \
      -n ${variant_set}
  
  ## Convert to plink format for faster computation with xerxes.
  ${params.poseidon_exec_dir}/trident genoconvert \
      --outFormat PLINK \
      -d .
  """
}