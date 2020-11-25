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

      --outdir [str]          The desired directory within which all output files will be placed. Subdirectories `data/<chrom_length>/<four_mN>` and `results/<chrom_length>/<four_mN>` will be created within this directory.

      --four_mN [float]           The scaled migration rate between non-diagonal neighbour populations used in the simulation.

      --chrom_length [float]  The length of the simulated chromosomes.

      --n_ind_per_pop [int]   The sample size of each of the 9 populations. All populations have an equal population and sample size.

      --max_ras_ac [int]      The maximum allele count to be considered for rare allele sharing.
      
      --knn [int]             The number of nearest neighbours to be used for KNN classification.
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
  publishDir "${params.outdir}/data/${params.chrom_length}/${params.four_mN}", mode: 'copy'
  memory '8GB'
  
  input:
  val chrom_name from Channel.of(1..20)
  
  output:
  tuple val(chrom_name), val("all"), path("all_vars_m${params.four_mN}_chr${chrom_name}.geno"), path("all_vars_m${params.four_mN}_chr${chrom_name}.snp"), path("all_vars_m${params.four_mN}_chr${chrom_name}.ind"), path("all_vars_m${params.four_mN}_chr${chrom_name}.indEach") into (ch_all_vars_datasets)
  tuple val(chrom_name), val("common"), path("common_vars_m${params.four_mN}_chr${chrom_name}.geno"), path("common_vars_m${params.four_mN}_chr${chrom_name}.snp"), path("common_vars_m${params.four_mN}_chr${chrom_name}.ind"), path("common_vars_m${params.four_mN}_chr${chrom_name}.indEach") into (ch_common_vars_datasets, ch_for_1240k_input_geno, ch_for_1240k_input_snp, ch_for_1240k_input_ind, ch_for_1240k_input_indEach)
  tuple val(chrom_name), val("rare"), path("all_vars_m${params.four_mN}_chr${chrom_name}.freqsum.gz") into ch_freqsum_dataset
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
  publishDir "${params.outdir}/data/${params.chrom_length}/${params.four_mN}", mode: 'copy'
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

process make_poplists {

  tag "m${params.four_mN}_l${params.chrom_length}"
  publishDir "${params.outdir}/results/${params.chrom_length}/${params.four_mN}/f3/poplists", mode: 'copy'
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

/*
  Concatenate channels of eigenstrat inputs by input type, then transpose to make into tables with 1 chromosome per row.
  Map to indicate the geno, ind, snp file.
  Then combine with poplists to make all possible pairs of inputs and poplists.
  These pairs will be used to submit f3 jobs for each chunk and variant type.
*/

ch_f3_input = ch_all_vars_datasets
                    .mix(ch_common_vars_datasets, ch_twelve_forty_vars_datasets)
                    .spread(ch_pairwise_poplists_f3)

process f3 {
//  conda 'bioconda::admixtools=6.0' // Added directly to environment.yml.
  tag "${variant_set}_chr${chrom_name}_m${params.four_mN}_l${params.chrom_length}"
  publishDir "${params.outdir}/results/${params.chrom_length}/${params.four_mN}/f3", mode: 'copy'
  memory '8GB'

  input:
  tuple chrom_name, variant_set, path(geno), path(snp), path(ind), path(indEach), path(poplist) from ch_f3_input

  output:
  tuple chrom_name, variant_set, path('*.out') into ch_f3_output

  script:
  """
  lineNr=\$(basename ${poplist} .txt | rev | cut -f1 -d '.' | rev)
  ## Make Parfile
  (echo -e "genotypename:\t${geno}"
  echo -e "snpname:\t${snp}"
  echo -e "indivname:\t${indEach}"
  echo -e "popfilename:\t${poplist}"
  echo -e "outgroupmode:\tYES"
  ) >parfile_f3.txt
  
  qp3Pop -p parfile_f3.txt >\$(basename ${geno} .geno).\$lineNr.out
  """
}

/* Group the f3_output channel by snp_set, so that one job is ran per variant set, and it contains all the output files associated with that variant set."*/
ch_f3_output
        .groupTuple(by: 1)
        .set{ ch_for_f3_matrix_generation }


process compile_F3_matrix {
  tag "${variant_set}_f3_matrix"
  publishDir "${params.outdir}/results/${params.chrom_length}/${params.four_mN}/similarity_matrices", mode: 'copy'
  memory '8GB'

  input:
  tuple chrom_name, variant_set, path(f3_logs) from ch_for_f3_matrix_generation

  output:
  tuple variant_set,path("${variant_set}_similarity_matrix.txt") into ch_eigenstrat_similarity_matrices

  script:
  """
  ${baseDir}/f3_to_similarity_matrix.py ${params.n_ind_per_pop} ${variant_set}_similarity_matrix.txt ${f3_logs}
  """
}

process run_Rascal {
  tag "m${params.four_mN}_chr${chrom_name}_l${params.chrom_length}"
  publishDir "${params.outdir}/results/${params.chrom_length}/${params.four_mN}/ras", mode: 'copy'
  memory '8GB'
  cpus 2

  input:
  tuple chrom_name, variant_set, path(freqsum_input) from ch_freqsum_dataset

  output:
  tuple chrom_name, variant_set, path("*.out") into ch_ras_output

  script:
  // First individual is 'ind0' so maximum individual number is 1 less than the product.
  num_inds= params.n_ind_per_pop * 9 - 1
  """
  ## Create list of all individuals from the number of individuals per population times 9.
  for ind1 in `seq 0 ${num_inds}`; do
    ## This will create a comma separated list with all individual names that has a trailing comma.
    ras_ind_list+="ind\${ind1},"
  done
  
  ## Ascertain in, and calculate RAS for all non-Ref populations (-a, -L, -c)
  zcat ${freqsum_input}  | ${params.rastools_dir}/RASCalculator.py --skipJackknife --details \
    -O ${variant_set}_m${params.four_mN}_chr${chrom_name}.out \
    -o Ref \
    -M ${params.max_ras_ac} \
    -c \${ras_ind_list%,} \
    -L \${ras_ind_list%,} \
    -a \${ras_ind_list%,} \
    -C 20
  """
}

/* Group the ras_output channel by snp_set, so that one job is ran per variant set, and it contains all the output files associated with that variant set."*/
ch_ras_output
        .groupTuple(by: 1)
        .set{ ch_for_ras_matrix_generation }

process compile_ras_matrix {
  tag "m${params.four_mN}_l${params.chrom_length}"
  publishDir "${params.outdir}/results/${params.chrom_length}/${params.four_mN}/similarity_matrices", mode: 'copy'
  memory '8GB'

  input:
  tuple chrom_name, variant_set, path(ras_logs) from ch_for_ras_matrix_generation
  
  output:
  tuple variant_set, path("rare_similarity_matrix.ac*.txt") into ch_rare_similarity_matrices
  tuple val("all.rascal"), path("all.rascal_similarity_matrix.txt") into ch_all_rascal_similarity_matrix
  
  script:
  """
  ${baseDir}/ras_to_similarity_matrices.py ${params.max_ras_ac} ${ras_logs}
  """
}


ch_eigenstrat_similarity_matrices
        .mix(ch_rare_similarity_matrices)
        .mix(ch_all_rascal_similarity_matrix)
/*        .dump(tag:"mixed similarity matrices")*/
        .into { ch_similarity_matrices; ch_similarity_matrices_for_Heatmap; ch_similarity_matrices_for_KNN}

process convert_to_distance_matrix{
  tag "m${params.four_mN}_l${params.chrom_length}"
  publishDir "${params.outdir}/results/${params.chrom_length}/${params.four_mN}/distance_matrices", mode: 'copy'
  memory '8GB'

  input:
    tuple variant_set, path(similarity_matrices) from ch_similarity_matrices

  output:
    tuple variant_set, path('*_distance_matrix*txt') into (ch_distance_matrices, ch_distance_matrices_for_KNN)

  script:
  """
  ${baseDir}/similarity_to_distance.py ${variant_set} ${similarity_matrices}
  """
}

// Filter similarity and distance matrices to exclude all rare variant matrices except the one with the maximum AC.
ch_distance_matrices
        .map{
          def snp_set = it[0]
          def files = it[1]
          def matrix = snp_set == "rare" ? files.find {it =~ /[A-Za-z0-9_]*.ac${params.max_ras_ac}.txt/ } : files
          
          [ matrix ]
        }
        .collect()
        .set { ch_for_MDS }

ch_similarity_matrices_for_Heatmap
        .map{
          def snp_set = it[0]
          def files = it[1]
          def matrix = snp_set == "rare" ? files.find {it =~ /[A-Za-z0-9_]*.ac${params.max_ras_ac}.txt/ } : files
  
          [ matrix ]
        }
        .collect()
        .set { ch_for_Heatmap }

process do_MDS {
  tag "m${params.four_mN}_l${params.chrom_length}"
  publishDir "${params.outdir}/plots/${params.chrom_length}/MDS/", mode: 'copy'
  memory '8GB'
  
  input:
  path(distance_matrices) from ch_for_MDS

  output:
  path('*pdf')

  script:
  """
  ${baseDir}/distance_to_MDS_plot.py ${params.four_mN} ${params.n_ind_per_pop} ${distance_matrices}
  """
}

process do_Heatmap {
  tag "m${params.four_mN}_l${params.chrom_length}"
  publishDir "${params.outdir}/plots/${params.chrom_length}/Heatmaps/", mode: 'copy'
  memory '8GB'
//  time '10m'
  
  input:
  path(similarity_matrices) from ch_for_Heatmap

  output:
  path('*pdf')

  script:
  """
  ${baseDir}/similarity_to_heatmap.py ${params.four_mN} ${params.n_ind_per_pop} ${similarity_matrices}
  """
}

process do_KNN {
  tag "m${params.four_mN}_${snp_set}_l${params.chrom_length}"
  publishDir "${params.outdir}/results/${params.chrom_length}/${params.four_mN}/KNN_classification", mode: 'copy'
  memory '8GB'
//  time '10m'
  
  input:
  tuple variant_set, path(similarity_matrices) from ch_distance_matrices_for_KNN

  output:
  path ('KNN*.txt') 
  path ('') into ch_dummy_plotting_delay // Output goes into a dummy channel that is only used to delay the plotting of all the KNN results across all four_mN runs thus far.

  script:
  """
  ${baseDir}/distance_to_KNN.py ${params.four_mN} ${params.n_ind_per_pop} ${variant_set} ${params.knn} ${similarity_matrices}
  """
}

/* Create a channel that picks up all KNN results from the same chrom_length regardless of four_mN value. 
  In that channel, mix the dummy delay channel and fiter for unique files to avoid any duplications. */
ch_for_KNN_plotting=Channel.fromPath("${params.outdir}/results/${params.chrom_length}/*/KNN_classification/KNN_K${params.knn}*.txt")
          .mix(ch_dummy_plotting_delay)
          .unique()
/*          .dump(tag:"for KNN Plot")*/

process plot_KNN {
  tag "l${params.chrom_length}"
  publishDir "${params.outdir}/plots/${params.chrom_length}/KNN_classification/", mode: 'copy'
  memory '8GB'
//  time '10m'

  input:
  path(knn_files) from ch_for_KNN_plotting.collect().dump()

  output:
  path('KNN_summary_plot*pdf')

  script:
  """
  ${baseDir}/plot_KNN.py ${params.chrom_length} ${knn_files}
  """
}