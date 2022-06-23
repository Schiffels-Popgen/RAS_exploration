#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  =========================================
  rare_variant_data_simulation_analysis v${workflow.manifest.version}
  =========================================
  Usage:

  The typical command for running the pipeline on sdag is as follows:

  nextflow run simulate_data.nf -profile cobra,mpcdf,conda,var_sim --four_mN 1 --chrom_length 1e6 --n_ind_per_pop 20 --max_ras_ac 5 --knn 5

  Mandatory arguments:
      -profile [str]          Institution or personal hardware config to use (e.g. standard, docker, singularity, conda, aws). Ask your system admin if unsure, or check documentation.

      --chrom_length [float]  The length of the simulated chromosomes. [100000000]

      --num_chroms [int]      The number of chromosomes to simulate. [20]

      --n_ind_per_pop [int]   The sample size of each of the 9 populations. All populations have an equal population and sample size. [500]
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
simulation_input = Channel.of(1..params.num_chroms)
                    .combine(Channel.fromList([1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0]))

def chrom_tag = "${params.chrom_length}x${params.num_chroms}"
process msprime{
  
  tag "n${params.n_ind_per_pop}_m${four_mN}_chr${chrom_name}_l${chrom_tag}"
  // publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}", mode: 'copy'
  memory '10GB'
  // executor 'local'

  input:
  tuple chrom_name, four_mN from simulation_input
  
  output:
  tuple val(four_mN), path("all_vars_m${four_mN}_chr${chrom_name}.geno") into ch_all_vars_geno
  tuple val(four_mN), path("all_vars_m${four_mN}_chr${chrom_name}.snp") into ch_all_vars_snp
  tuple val(four_mN), path("all_vars_m${four_mN}_chr${chrom_name}.ind") into ch_all_vars_ind
  tuple val(four_mN), path("common_vars_m${four_mN}_chr${chrom_name}.geno") into (ch_common_vars_geno, ch_common_vars_geno_for_1240k )
  tuple val(four_mN), path("common_vars_m${four_mN}_chr${chrom_name}.snp") into (ch_common_vars_snp, ch_common_vars_snp_for_1240k )
  tuple val(four_mN), path("common_vars_m${four_mN}_chr${chrom_name}.ind") into (ch_common_vars_ind, ch_common_vars_ind_for_1240k )
  path("variant_counts.m${four_mN}_chr${chrom_name}.txt") 
  
  script:
  """
  ${baseDir}/simulation.py -c ${params.chrom_length} -n ${chrom_name} -s ${params.n_ind_per_pop} -m ${four_mN} #-o /projects1/MICROSCOPE/rarevar_sim_study/data/
  ${baseDir}/add_ref_and_sort_eigenstrat.sh ${chrom_name} ${four_mN} .
  """
}

process make_1240k{
  tag "n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
  publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}", mode: 'copy'
  memory '1GB'
  // executor 'local'

  input:
  tuple val(four_mN), path(genos) from ch_common_vars_geno_for_1240k.groupTuple(by:0)//.dump(tag:"make_1240k_input_geno")
  tuple val(four_mN), path(snps) from ch_common_vars_snp_for_1240k.groupTuple(by:0)//.dump(tag:"make_1240k_input_snp")
  tuple val(four_mN), path(inds) from ch_common_vars_ind_for_1240k.groupTuple(by:0)//.dump(tag:"make_1240k_input_ind")

  output:
  tuple val("1:${params.num_chroms}"), val(four_mN), val("twelve_forty_vars"),path("twelve_forty_vars_m${four_mN}_l${params.chrom_length}.geno"),path("twelve_forty_vars_m${four_mN}_l${params.chrom_length}.snp"),path("twelve_forty_vars_m${four_mN}_l${params.chrom_length}.ind") into ch_twelve_forty_vars_datasets

  script:
  """
  ## First concatenate the snps and genos by chromosome
  for chrom_name in {1..${params.num_chroms}}; do
    cat common_vars_m${four_mN}_chr\${chrom_name}.geno  >> concat_genos
    cat common_vars_m${four_mN}_chr\${chrom_name}.snp  >> concat_snps
  done
  
  ## Then paste side by side, use shuf to randomly subsample 1200K sites, sort by original line numbers, and remove leading spaces (added by cat -n)
  paste -d " " concat_snps concat_genos | cat -n | shuf -n 1200000 | sort -nk1 | sed -e 's/^[ ]*//' >banana.txt
  
  ## Extract the .snp and .geno file from the subsample. 
  cut -f 1 -d " " banana.txt | cut -f 2- >twelve_forty_vars_m${four_mN}_l${params.chrom_length}.snp
  cut -f 2 -d " " banana.txt >twelve_forty_vars_m${four_mN}_l${params.chrom_length}.geno
  
  ## Copy the original ind files 
  cp ${inds[0]} twelve_forty_vars_m${four_mN}_l${params.chrom_length}.ind
  """
}

process merge_chroms_common_vars {
  tag "n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
  publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/", mode: 'copy'
  memory '50MB'
  cpus 1
  // executor 'local'

  input:
  tuple val(four_mN), path(genos) from ch_common_vars_geno.groupTuple(by:0)//.dump(tag:"merge_common_input_geno")
  tuple val(four_mN), path(snps) from ch_common_vars_snp.groupTuple(by:0)//.dump(tag:"merge_common__input_snp")
  tuple val(four_mN), path(inds) from ch_common_vars_ind.groupTuple(by:0)//.dump(tag:"merge_common__input_ind")

  output:
  tuple val("1:${params.num_chroms}"), val(four_mN), val("common_vars"),path("common_vars_m${four_mN}_l${params.chrom_length}.geno"),path("common_vars_m${four_mN}_l${params.chrom_length}.snp"),path("common_vars_m${four_mN}_l${params.chrom_length}.ind") into ch_common_vars_datasets

  script:
  """
  ## First concatenate the snps and genos by chromosome
  for chrom_name in {1..${params.num_chroms}}; do
    cat common_vars_m${four_mN}_chr\${chrom_name}.geno  >> common_vars_m${four_mN}_l${params.chrom_length}.geno
    cat common_vars_m${four_mN}_chr\${chrom_name}.snp  >> common_vars_m${four_mN}_l${params.chrom_length}.snp
  done

  ## Copy the original ind files 
  cp ${inds[0]} common_vars_m${four_mN}_l${params.chrom_length}.ind
  """
}


process merge_chroms_all_vars {
  tag "n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
  publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/", mode: 'copy'
  memory '50MB'
  cpus 1
  // executor 'local'

  input:
  tuple val(four_mN), path(genos) from ch_all_vars_geno.groupTuple(by:0)//.dump(tag:"merge_all_input_geno")
  tuple val(four_mN), path(snps) from ch_all_vars_snp.groupTuple(by:0)//.dump(tag:"merge_all_input_snp")
  tuple val(four_mN), path(inds) from ch_all_vars_ind.groupTuple(by:0)//.dump(tag:"merge_all_input_ind")

  output:
  tuple val("1:${params.num_chroms}"), val(four_mN), val("all_vars"),path("all_vars_m${four_mN}_l${params.chrom_length}.geno"),path("all_vars_m${four_mN}_l${params.chrom_length}.snp"),path("all_vars_m${four_mN}_l${params.chrom_length}.ind") into ch_all_vars_datasets

  script:
  """
  ## First concatenate the snps and genos by chromosome
  for chrom_name in {1..${params.num_chroms}}; do
    cat all_vars_m${four_mN}_chr\${chrom_name}.geno  >> all_vars_m${four_mN}_l${params.chrom_length}.geno
    cat all_vars_m${four_mN}_chr\${chrom_name}.snp  >> all_vars_m${four_mN}_l${params.chrom_length}.snp
  done

  ## Copy the original ind files 
  cp ${inds[0]} all_vars_m${four_mN}_l${params.chrom_length}.ind
  """
}

ch_all_vars_datasets
    .mix(ch_common_vars_datasets, ch_twelve_forty_vars_datasets)
    // .dump(tag:"datasets")
    .into { ch_datasets; ch_datasets_for_subsampling}

process create_poseidon_packages {
  tag "${variant_set}_n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
  publishDir "${baseDir}/../data/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/poseidon", mode: 'copy'
  memory '1GB'
  cpus 1
  // executor 'local'

  input:
  tuple chroms, four_mN, variant_set, path(geno), path(snp), path(ind) from ch_datasets

  output:
  tuple path("${variant_set}/${pkg_name}.bed"), path("${variant_set}/${pkg_name}.bim"), path("${variant_set}/${pkg_name}.fam"), path("${pkg_name}.geno"), path("${pkg_name}.snp"), path("${pkg_name}.ind"), path("${variant_set}/POSEIDON.yml")

  script:
  pkg_name="${variant_set}_m${four_mN}_l${params.chrom_length}"
  """
  ## trident creates the package within the work directory, and nextflow is responsible for putting it in the data dir.
  ## Skip genoconvert by using forge to covert to plink on the fly (v1.0.0.0+). 
  ##    The original data is then copied over by nextflow into the package dir to have also eigenstrat format.

    ${params.poseidon_exec_dir}/trident forge -d . \
        -p ${pkg_name}.geno \
        -o ${variant_set} \
        -n ${pkg_name} \
        --minimal ## No janno or bib file
  """
}

ch_subsampling_dataset = Channel.fromList( [10,20,50,100,200] )
          .filter { it < params.n_ind_per_pop } // When n_ind_per_pop is below 500, only subsample to lower values than the set one.
          .combine(ch_datasets_for_subsampling)

process subsample_datasets{
  tag "${variant_set}_m${four_mN}_n${n}_l${chrom_tag}"
  publishDir "${baseDir}/../data/n${n}/${chrom_tag}/${four_mN}/poseidon", mode: 'copy'
  memory '1GB'
  cpus 1
  // executor 'local'

  input:
  tuple val(n), chroms, four_mN, variant_set, path(geno), path(snp), path(ind) from ch_subsampling_dataset

  output:
  tuple path("${variant_set}/${pkg_name}.bed"), path("${variant_set}/${pkg_name}.bim"), path("${variant_set}/${pkg_name}.fam"), path("downsampled_data/${pkg_name}.geno"), path("downsampled_data/${pkg_name}.snp"), path("downsampled_data/${pkg_name}.ind"), path("${variant_set}/POSEIDON.yml")
  
  script:
  pkg_name="${variant_set}_m${four_mN}_l${params.chrom_length}"
  """
  mkdir -p downsampled_data

  ## Make list of coordinates needed for cutting genos out
  gidxs=''  ## Genotype indexes
  iidxs=''  ## Individual indexes
  let start=1
  for p in {0..8}; do
    let end=\$(( \${start} + ${n} - 1 ))
    gidxs+="\${start}-\${end},"
    iidxs+="\${start},\${end}p; "
    let start+=${params.n_ind_per_pop}
  done
  ## include "Ref" individual
  let ref=\$(( ${params.n_ind_per_pop} * 9 + 1 )) ## It is always added at the end of the ind file
  gidxs+="\${ref}"
  iidxs+="\${ref}p;"
  
  ## Copy .snp intact
  cp ${snp} downsampled_data/
  ## Downsample .geno
  cut -c \${gidxs} ${geno} >downsampled_data/${pkg_name}.geno
  ## Dowsample .ind
  sed -n "\${iidxs}" ${ind} >downsampled_data/${pkg_name}.ind

  ## Skip genoconvert by using forge to covert to plink on the fly (v1.0.0.0+). 
  ##    The original data is then copied over by nextflow into the package dir to have also eigenstrat format.
  ${params.poseidon_exec_dir}/trident forge -d . \
      -p downsampled_data/${pkg_name}.geno \
      -o ${variant_set} \
      -n ${pkg_name} \
      --minimal ## No janno or bib file
  """
}
