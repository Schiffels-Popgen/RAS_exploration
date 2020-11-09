#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  =========================================
  rare_variant_simulation_analysis v${workflow.manifest.version}
  =========================================
  Usage:

  The typical command for running the pipeline on sdag is as follows:

  nextflow run pipeline.nf -profile conda,shh,sdag --4mN 1 --chrom_length 1e6 --n_ind_per_pop 20

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
  publishDir "/projects1/MICROSCOPE/rarevar_sim_study/data/${params.chrom_length}/${params.four_mN}", mode: 'copy'
  queue "long"
  memory '8GB'
  
  input:
  val chrom_name from Channel.of(1..20)
  
  output:
  path("all_vars_m${params.four_mN}_chr${chrom_name}.geno") into ch_all_vars_geno_for_f3
  path("all_vars_m${params.four_mN}_chr${chrom_name}.snp") into ch_all_vars_snp_for_f3
  path("all_vars_m${params.four_mN}_chr${chrom_name}.ind") into ch_all_vars_ind_for_f3
  path("all_vars_m${params.four_mN}_chr${chrom_name}.indEach") into ch_all_vars_indEach_for_f3
  path("common_vars_m${params.four_mN}_chr${chrom_name}.geno") into (ch_common_vars_geno_for_f3, ch_common_vars_geno_for_1240k )
  path("common_vars_m${params.four_mN}_chr${chrom_name}.snp") into (ch_common_vars_snp_for_f3, ch_common_vars_snp_for_1240k )
  path("common_vars_m${params.four_mN}_chr${chrom_name}.ind") into (ch_common_vars_ind_for_f3, ch_common_vars_ind_for_1240k )
  path("common_vars_m${params.four_mN}_chr${chrom_name}.indEach") into (ch_common_vars_indEach_for_f3, ch_common_vars_indEach_for_1240k )
  path("variant_counts.m${params.four_mN}_chr${chrom_name}.txt") 
  
  script:
  """
  /projects1/MICROSCOPE/rarevar_sim_study/scripts/simulation.py -c ${params.chrom_length} -n ${chrom_name} -s ${params.n_ind_per_pop} -m ${params.four_mN} #-o /projects1/MICROSCOPE/rarevar_sim_study/data/
  /projects1/MICROSCOPE/rarevar_sim_study/scripts/add_ref_and_sort_eigenstrat.sh ${chrom_name} ${params.four_mN} .
  """
}


process make_1240k{

  tag "m${params.four_mN}_chr${chrom_name}_l${params.chrom_length}"
  publishDir "/projects1/MICROSCOPE/rarevar_sim_study/data/${params.chrom_length}/${params.four_mN}", mode: 'copy'
  queue "short"
  memory '1GB'

  input:
  path genos from ch_common_vars_geno_for_1240k.collect()
  path snps from ch_common_vars_snp_for_1240k.collect()
  path inds from ch_common_vars_ind_for_1240k.collect()
  path indEachs from ch_common_vars_indEach_for_1240k.collect()
  
  output:
  path("twelve_forty_vars_m${params.four_mN}.geno") into ch_twelve_forty_vars_geno_for_f3
  path("twelve_forty_vars_m${params.four_mN}.snp") into ch_twelve_forty_vars_snp_for_f3
  path("twelve_forty_vars_m${params.four_mN}.ind") into ch_twelve_forty_vars_ind_for_f3
  path("twelve_forty_vars_m${params.four_mN}.indEach") into ch_twelve_forty_vars_indEach_for_f3

  script:
  """
  ## First concatenate the snps and genos by chromosome
  for chrom_name in {1..20}; do
    cat common_vars_m${params.four_mN}_chr\${chrom_name}.geno  >> concat_genos
    cat common_vars_m${params.four_mN}_chr\${chrom_name}.snp  >> concat_snps
  done
  
  ## Then paste side by side, use shuf to randomly subsample 1200K sites
  paste -d " " concat_snps concat_genos | cat -n | shuf -n 1200000 | sort -nk1 >banana.txt
  
  ## Extract the .snp and .geno file from the subsample. 
  cut -f 2-7 -d " " banana.txt >twelve_forty_vars_m${params.four_mN}.snp
  cut -f 8 -d " " banana.txt >twelve_forty_vars_m${params.four_mN}.geno
  
  ## Copy the original ind files 
  cp ${inds[0]} twelve_forty_vars_m${params.four_mN}.ind
  cp ${indEachs[0]} twelve_forty_vars_m${params.four_mN}.indEach
  """

}

process make_poplists {

  tag "m${params.four_mN}_l${params.chrom_length}"
  publishDir "/projects1/MICROSCOPE/rarevar_sim_study/results/f3/poplists", mode: 'copy'
  queue "short"
  memory '1GB'

  output:
  path("pairwise_poplist.*.txt") into ch_pairwise_poplists_f3
  file("pairwise_poplist.txt")

  script:
  num_inds = params.n_ind_per_pop * 9 - 1

  // 1500 lines per job.
  lines_per_job = 1500
  """
  ## Create overall poplist
  for i in \$(seq 0 1 ${num_inds}); do
    lower_bound=\$(expr \${i} + 1)
    for j in \$(seq \${lower_bound} 1 ${num_inds}); do
      echo -e "ind\${i}\tind\${j}\tRef"
    done
  done >pairwise_poplist.txt

  ## Split into chunks for paralellisation
  let all_lines=\$(wc -l pairwise_poplist.txt | cut -d " " -f 1)
  let start_line=1
  while [[ \${start_line} -lt \${all_lines} ]]; do
    let end_line=\$(expr \${start_line} + ${lines_per_job})
    sed -n -e \${start_line},\${end_line}p pairwise_poplist.txt > pairwise_poplist.\${start_line}.txt
    let start_line=\$(expr \${end_line} + 1)
  done
  """
}

/*process f3 {
  tag "m${params.four_mN}_l${params.chrom_length}"
  publishDir "/projects1/MICROSCOPE/rarevar_sim_study/results/f3/poplists", mode: 'copy'
  queue "short"
  memory '8GB'

  input:
    path("pairwise_poplist.*") into ch_pairwise_poplists_f3.collect()


  output:

}*/
/*ch_all_vars_for_tweaking.view()*/
