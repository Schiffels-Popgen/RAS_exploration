#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  =========================================
  rare_variant_simulation_analysis v${workflow.manifest.version}
  =========================================
  Usage:

  The typical command for running the pipeline on sdag is as follows:

  nextflow run analyse_data.nf -profile cobra,mpcdf,conda,var_sim --four_mN 1 --chrom_length 1e6 --n_ind_per_pop 20 --max_ras_ac 5 --knn 5

  Mandatory arguments:
      -profile [str]          Institution or personal hardware config to use (e.g. standard, docker, singularity, conda, aws). Ask your system admin if unsure, or check documentation.

      --four_mN [float]           The scaled migration rate between non-diagonal neighbour populations used in the simulation.

      --chrom_length [float]  The length of the simulated chromosomes. [100000000]

      --num_chroms [int]      The number of chromosomes to simulate. [20]

      --total_inds_per_pop [int]  The number of individuals in the original large dataset, prior to subsampling. [500]

      --n_ind_per_pop [int]   The sample size of each of the 9 populations. All populations have an equal population and sample size.

      --max_ras_ac [int]      The maximum allele count to be considered for rare allele sharing.
      
      --knn [int]             The number of nearest neighbours to be used for KNN classification.
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

// Parameters to refer to within multiple processes. 
def chrom_tag = "${params.chrom_length}x${params.num_chroms}"

ch_input_datasets = Channel.fromList(["all_vars", "common_vars", "twelve_forty_vars", "rare_vars"])
      .combine(Channel.fromList( [1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0,500.0]))
      // .combine(Channel.fromList( [10, 20, 50, 100, 200, 500] ))
      .map{
        def variant_set = it[0]
        def pkg_name = "all_vars" // Use 'all_vars' package for all runs, and rely on xerxes to do the allele frequency ascertainment in each case
        // def pkg_name = variant_set == "rare_vars" ? "all_vars" : variant_set
        def four_mN = it[1]
        // def n_ind_per_pop = it[2]
        def package_dir = "${baseDir}/../data/n${params.n_ind_per_pop}/${chrom_tag}/"+it[1]+"/poseidon/"+pkg_name
        // Pick up the actual files to ensure correct resume behaviour
        def bed = "${package_dir}/${pkg_name}_m${four_mN}_l${params.chrom_length}.bed"
        def bim = "${package_dir}/${pkg_name}_m${four_mN}_l${params.chrom_length}.bim"
        def fam = "${package_dir}/${pkg_name}_m${four_mN}_l${params.chrom_length}.fam"

        // [ variant_set, four_mN, n_ind_per_pop, package_dir, bed, bim, fam ]
        [ variant_set, four_mN, package_dir, bed, bim, fam ]
      }
      .branch{
          all: it[0] == "all_vars"
          other: true
      }

// Set of minimum(top) and maximum(bottom) AFs to run on all_vars dataset. 0.0-1.0 included also as final set.
ch_af_cutoffs = Channel.from( [0.01, 0.02, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.8, 0.0] )
          .merge(Channel.from( [0.02, 0.05, 0.10, 0.15, 0.2 , 0.3, 0.4, 0.5, 0.8, 1.0, 1.0]))

ch_input_datasets.other.map{
    variant_set, four_mN, package_dir, bed, bim, fam ->
    def minAF = variant_set == "rare_vars" ? 2 : 0.05
    def maxAF = variant_set == "rare_vars" ? params.max_ras_ac : 0.95
    [ variant_set, four_mN, package_dir, bed, bim, fam , minAF, maxAF ]
  }
  .into{ ch_non_all_vars_for_pairwise_ras_prep; ch_non_all_vars_for_ras_prep }

ch_input_datasets.all.combine(ch_af_cutoffs)
      .into{ch_all_vars_for_pairwise_ras_prep; ch_all_vars_for_ras_prep}

// Mix all and non-all vars with just 0-1 AF bin for all_vars
ch_all_vars_for_pairwise_ras_prep
  .filter( {it[6] == 0.0 && it[7] == 1.0} )
  .mix( ch_non_all_vars_for_pairwise_ras_prep )
  .dump(tag:"input_pairwise_ras")
  .into{ ch_input_xerxes_pairwise_ras }

// Mix all and non-all vars with many AF bins
ch_all_vars_for_ras_prep.mix( ch_non_all_vars_for_ras_prep )
  .dump(tag:"input_ras")
  .into { ch_input_for_xerxes_ras }

// // This might be overkill, but split f3/f4 runs to only run statistics with 1 individual at position A. Then merge all logfiles together later
// // First create a channel of the individuals in this (potentially subsampled) dataset.
// ch_individuals = Channel.from( 0..8 )
//                       .combine(Channel.from( 0..<params.n_ind_per_pop )) //< means the actual value of params.n_ind_per_pop is not in the range.
//                       .map{
//                               first_ind = it[0] * params.total_inds_per_pop
//                               actual_ind= it[1] + first_ind
//                               actual_ind
//                           }

// ch_for_combining
//   .combine(ch_individuals) 
//   .into{ ch_input_for_xerxes_f3; ch_input_for_xerxes_f4; }

process xerxes_pairwise_ras {
  tag "${variant_set}_n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}_m${minAF}_M${maxAF}"
  publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/xerxes_ras", mode: 'copy'
  memory '64GB'
  cpus 1
//  executor 'local'

  input:
  tuple val(variant_set), val(four_mN), val(package_dir), path(bed), path(bim), path(fam), val(minAF), val(maxAF) from ch_input_xerxes_pairwise_ras
  // tuple val(variant_set), val(four_mN), val(n_ind_per_pop), val(package_dir), path(bed), path(bim), path(fam) from ch_input_xerxes_pairwise_ras

  output:
  tuple variant_set, four_mN, path("pairwise_ras_table_${variant_set}_m${minAF}_M${maxAF}.out") into (ch_xerxes_pairwise_ras_output_for_matrix)
  file "pairwise_popConfigFile_${variant_set}.txt"
  // file "pairwise_blockTableFile_${variant_set}_m${minAF}_M${maxAF}.txt"

  script:
  def freq_cutoffs = variant_set == "rare_vars" ? "--minAC ${minAF} --maxAC ${maxAF}" : minAF == 0.0 && maxAF == 1.0 ? "--noMinFreq --noMaxFreq" : "--minFreq ${minAF} --maxFreq ${maxAF}"
  def partial_run = variant_set == "twelve_forty_vars" ? "--maxSnps 12000000" : ''
  def inds_to_run = params.n_ind_per_pop <= 20 ? params.n_ind_per_pop : 20 // Use 20 individuals per population for pairwise matrix. if fewer than 10 are there, use all.
  """
  ## Create popConfigFile. Each pop consists of all but the first two individuals. Pairwise matrix only calculated for 20 individuals
  echo "groupDefs:" >pairwise_popConfigFile_${variant_set}.txt

  ## Define rights and lefts
  echo "popLefts:" >>pairwise_popConfigFile_${variant_set}.txt
  for p in \$(seq 0 1 8); do
    idx=\$(( \${p} * ${params.total_inds_per_pop}))
    for j in \$(seq \${idx} 1 \$(( \${idx} + ${inds_to_run} -1)) ); do
      echo "  - <ind\${j}>" >>pairwise_popConfigFile_${variant_set}.txt
    done
  done

  echo "popRights:" >>pairwise_popConfigFile_${variant_set}.txt
  for p in \$(seq 0 1 8); do
    idx=\$(( \${p} * ${params.total_inds_per_pop}))
    for j in \$(seq \${idx} 1 \$(( \${idx} + ${inds_to_run} -1)) ); do
      echo "  - <ind\${j}>" >>pairwise_popConfigFile_${variant_set}.txt
    done
  done

  echo "outgroup: <Ref>" >>pairwise_popConfigFile_${variant_set}.txt

  ## Run ras
  ${params.poseidon_exec_dir}/xerxes ras -d ${package_dir}/ \
    --popConfigFile pairwise_popConfigFile_${variant_set}.txt \
    ${freq_cutoffs} \
    ${partial_run} \
    -j CHR \
    > pairwise_ras_table_${variant_set}_m${minAF}_M${maxAF}.out
  """
}

process xerxes_ras {
  tag "${variant_set}_n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}_m${minAF}_M${maxAF}"
  publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/xerxes_ras", mode: 'copy'
  memory '8GB'
  cpus 1

  input:
  tuple val(variant_set), val(four_mN), val(package_dir), path(bed), path(bim), path(fam), val(minAF), val(maxAF) from ch_input_for_xerxes_ras

  output:
  tuple variant_set, four_mN, path("ras_table_${variant_set}_m${minAF}_M${maxAF}.out") into (ch_xerxes_ras_output_for_matrix, ch_xerxes_ras_for_rasta)
  file "popConfigFile_${variant_set}.txt"
  // file "blockTableFile_${variant_set}.txt"
  file "rasta_table_${variant_set}_m${minAF}_M${maxAF}.out"

  script:
  def freq_cutoffs = variant_set == "rare_vars" ? "--minAC ${minAF} --maxAC ${maxAF}" : minAF == 0.0 && maxAF == 1.0 ? "--noMinFreq --noMaxFreq" : "--minFreq ${minAF} --maxFreq ${maxAF}" 
  def partial_run = variant_set == "twelve_forty_vars" ? "--maxSnps 12000000" : ''

  """
  ## Create population definitions file. Each pop consists of all but the first two individuals. The first two individuals of a pop is used as a left pop.
  echo "groupDefs:" >popConfigFile_${variant_set}.txt
  leftpops=()
  rightpops=()
  for pop in {0..8}; do
    excluded_ind='-<ind'\$((${params.total_inds_per_pop} * \${pop}))'>'
    excluded_ind2='-<ind'\$((${params.total_inds_per_pop} * \${pop} + 1))'>'
    echo -e "  Pop\${pop}Rest: [ Pop\${pop},\${excluded_ind},\${excluded_ind2} ]" >>popConfigFile_${variant_set}.txt
    leftpops+=(\${excluded_ind#"-"}) ## Add excluded ind to lefts without the leading "-"
    leftpops+=(\${excluded_ind2#"-"}) ## Add excluded ind to lefts without the leading "-"
    rightpops+=(Pop\${pop}Rest) ## Add created population into rights
  done

  ## Define rights and lefts
  echo "popLefts:" >>popConfigFile_${variant_set}.txt
  for left in \${leftpops[@]}; do
    echo "  - \${left}" >>popConfigFile_${variant_set}.txt
  done

  echo "popRights:" >>popConfigFile_${variant_set}.txt
  for right in \${rightpops[@]}; do
    echo "  - \${right}" >>popConfigFile_${variant_set}.txt
  done

  echo "outgroup: <Ref>" >>popConfigFile_${variant_set}.txt

  ## Run ras
  ${params.poseidon_exec_dir}/xerxes ras -d ${package_dir}/ \
    --popConfigFile popConfigFile_${variant_set}.txt \
    ${freq_cutoffs} \
    ${partial_run} \
    -j CHR \
    --f4TableOutFile rasta_table_${variant_set}_m${minAF}_M${maxAF}.out \
    > ras_table_${variant_set}_m${minAF}_M${maxAF}.out
  """
}

// process xerxes_f4 {
//   tag "ind${ind_id} ${variant_set} n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
//   publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/xerxes_f3/individual_f3", mode: 'copy'
//   memory '8GB'
//   cpus 1

//   input:
//   tuple val(variant_set), val(four_mN), val(package_dir), path(bed), path(bim), path(fam), val(ind_id) from ch_input_for_xerxes_f3

//   output:
//   tuple variant_set, four_mN, ind_id, path("*.out") into (ch_xerxes_individual_f3)

//   script:
//   """
//   ## Function to join by multi-character sep
//   function join_by {
//     local d=\${1-} f=\${2-}
//     if shift 2; then
//       printf %s "\${f}" "\${@/#/\${d}}"
//     fi
//   }

//   ## Create population definitions file (requires xerxes 0.2.0.0)
//   leftpops=()
//   rightpops=()
//   ascertainOn=''
//   echo "groupDefs:" >>popConfigFile_${ind_id}.txt
//   for pop in {0..8}; do
//     excluded_ind='-<ind'\$((${params.total_inds_per_pop} * \${pop}))'>'
//     excluded_ind2='-<ind'\$((${params.total_inds_per_pop} * \${pop} + 1))'>'
//     echo -e "  Pop\${pop}Rest:[ \"Pop\${pop}\",\"\${excluded_ind}\",\"\${excluded_ind2}\" ]" >>popConfigFile_${ind_id}.txt
//     leftpops+=(\${excluded_ind#"-"}) ## Add excluded ind to lefts without the leading "-"
//     leftpops+=(\${excluded_ind2#"-"}) ## Add second excluded ind to lefts without the leading "-"
//     rightpops+=(Pop\${pop}Rest) ## Add created population into rights
//     ascertainOn+='"Pop\${pop}Rest", '
//   done
//   ## One more group definition of all pops for ascertainment
//   echo "AscertainOn: [ \${ascertainOn%, } ]" >>popConfigFile_${ind_id}.txt

//   echo "fstats:" >>popConfigFile_${ind_id}.txt
//   echo "- type: F4" >>popConfigFile_${ind_id}.txt
//   echo "  a: [ \" \$(join_by '\",\"' \${leftpops[@]}) \" ]" >>popConfigFile_${ind_id}.txt
//   echo "  b: [ \" \$(join_by '\",\"' \${rightpops[@]}) \" ]" >>popConfigFile_${ind_id}.txt
//   echo "  c: [ \"Ref\" ]" >>popConfigFile_${ind_id}.txt
//   echo "  ascertainment:" >>popConfigFile_${ind_id}.txt
//   echo "    outgroup: \"Ref\"" >>popConfigFile_${ind_id}.txt
//   echo "    reference: \"AscertainOn\"" >>popConfigFile_${ind_id}.txt
//   echo "    lower: 0.05 >>popConfigFile_${ind_id}.txt
//   echo "    upper: 0.95 >>popConfigFile_${ind_id}.txt



//   """
// }

// process xerxes_f3 {
//   tag "ind${ind_id} ${variant_set} n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
//   publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/xerxes_f3/individual_f3", mode: 'copy'
//   memory '8GB'
//   cpus 1
// //  executor 'local'

//   input:
//   tuple val(variant_set), val(four_mN), val(package_dir), path(bed), path(bim), path(fam), val(ind_id) from ch_input_for_xerxes_f3

//   output:
//   tuple variant_set, four_mN, ind_id, path("*.out") into (ch_xerxes_individual_f3)

//   script:
//   """
//   ## Create statFile
//   ## Create overall poplist
//   ## Because the individuals are downsampled from the complete simulation, they are not sequential between pops.
//   ##    For every pop, reseet index to next 500 multiple and then count the first X individuals from that.
//   for p in \$(seq 0 1 8); do
//     idx=\$(( \${p} * ${params.total_inds_per_pop}))
//     for j in \$(seq \${idx} 1 \$(( \${idx}+${params.n_ind_per_pop}-1)) ); do
//       echo "F3vanilla(<ind${ind_id}>, <ind\${j}>, Ref)"
//     done
//   done >pairwise_poplist_${variant_set}.ind${ind_id}.txt

//   ## Remove the last byte (i.e. the final newline character) from the statFile
//   truncate -s -1 pairwise_poplist_${variant_set}.ind${ind_id}.txt

//   ## Run f3
//   ${params.poseidon_exec_dir}/xerxes fstats -d ${package_dir}/ \
//     --statFile pairwise_poplist_${variant_set}.ind${ind_id}.txt \
//     -j CHR \
//     > f3_table_${variant_set}.ind${ind_id}.out
//   """
// }

// process xerxes_f4 {
//   tag "ind${ind_id} ${variant_set} n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
//   publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/xerxes_f3/individual_f4", mode: 'copy'
//   memory '16GB'
//   cpus 1
// //  executor 'local'

//   input:
//   tuple val(variant_set), val(four_mN), val(package_dir), path(bed), path(bim), path(fam), val(ind_id) from ch_input_for_xerxes_f4

//   output:
//   tuple variant_set, four_mN, ind_id, path("*.out") into (ch_xerxes_individual_f4)

//   script:
//   """
//   ## Create statFile
//   for b in \$(seq 0 1 8); do
//     for c in \$(seq 0 1 8); do
//       if [[ ! "\${b}" == "\${c}" ]]; then
//         echo -e "F4(<ind${ind_id}>, Pop\${b}, Pop\${c}, Ref)"
//       fi
//     done
//   done >pairwise_poplist_${variant_set}.ind${ind_id}.txt

//   ## Remove the last byte (i.e. the final newline character) from the statFile
//   truncate -s -1 pairwise_poplist_${variant_set}.ind${ind_id}.txt

//   ## Run f4
//   ${params.poseidon_exec_dir}/xerxes fstats -d ${package_dir}/ \
//     --statFile pairwise_poplist_${variant_set}.ind${ind_id}.txt \
//     -j CHR \
//     -f f4_table_${variant_set}.ind${ind_id}.out
//   """
// }

// ch_xerxes_f3_for_merging = ch_xerxes_individual_f3
//                             .groupTuple(by:[0,1]) // each element is now all individual f4 outputs across a variant_set/four_mN combination
//                             .dump(tag:"ch_xerxes_f3_for_merging")


// ch_xerxes_f4_for_merging = ch_xerxes_individual_f4
//                             .groupTuple(by:[0,1]) // each element is now all individual f4 outputs across a variant_set/four_mN combination
//                             .dump(tag:"ch_xerxes_f4_for_merging")

// process merge_f3 {
//   tag "${variant_set} n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
//   publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/xerxes_f3", mode: 'copy'
//   memory '1GB'
//   cpus 1
//   executor 'local' // No need to batch a simple cat command

//   input:
//   tuple variant_set, four_mN, ind_id, path(individual_f3_fn) from ch_xerxes_f3_for_merging

//   output:
//   tuple  variant_set, four_mN,  path("f3_table_${variant_set}.out") into (ch_xerxes_fstats_output_for_matrix, ch_xerxes_fstats_for_rasta)

//   script:
//   """
  
//   cat  f3_table_${variant_set}.*.out > f3_table_${variant_set}.out
//   """
// }

// process merge_f4 {
//   tag "${variant_set} n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
//   publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/xerxes_f3", mode: 'copy'
//   memory '1GB'
//   cpus 1
//   executor 'local' // No need to batch a simple cat command

//   input:
//   tuple variant_set, four_mN, ind_id, path(individual_f4_fn) from ch_xerxes_f4_for_merging

//   output:
//   tuple  variant_set, four_mN,  path("f4_table_${variant_set}.out") into (ch_xerxes_f4_output)

//   script:
//   """
//   cat  f4_table_${variant_set}.*.out > f4_table_${variant_set}.out
//   """
// }

// // Mix pairwise f3 and ras outputs
// ch_input_for_make_similarity_matrices = ch_xerxes_fstats_output_for_matrix
//     .mix(ch_xerxes_pairwise_ras_output_for_matrix)
//     .dump(tag:"ch_input_for_make_similarity_matrices")

// process make_similarity_matrix {
//   tag "${variant_set} n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
//   publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/similarity_matrices", mode: 'copy'
//   memory '1GB'
//   cpus 1
//   executor 'local'

//   input:
//   tuple variant_set, four_mN, file(pairwise_table) from ch_input_for_make_similarity_matrices

//   output:
//   tuple variant_set, four_mN, file("*_similarity_matrix.txt") into ch_similarity_matrices

//   script:
//   """
//   ${baseDir}/make_similarity_matrices.R ${variant_set} ${pairwise_table}
//   """
// }

// process make_distance_matrices {
//   tag "${variant_set} n${params.n_ind_per_pop}_m${four_mN}_l${chrom_tag}"
//   publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/distance_matrices", mode: 'copy'
//   memory '1GB'
//   cpus 1
//   executor 'local'

//   input:
//   tuple variant_set, four_mN, file(similarity_matrix) from ch_similarity_matrices

//   output:
//   tuple variant_set, four_mN, file("*_distance_matrix.txt") into (ch_distance_matrices_for_knn, ch_distance_matrices_for_mds)

//   script:
//   """
//   ${baseDir}/similarity_to_distance.py ${variant_set} ${similarity_matrix}
//   """
// }

// process do_MDS {
//   tag "n${params.n_ind_per_pop}_m${params.four_mN}_l${chrom_tag}"
//   publishDir "${baseDir}/../plots/n${params.n_ind_per_pop}/${chrom_tag}/MDS/", mode: 'copy'
//   memory '8GB'
//   cpus 1
//   executor 'local'
  
//   input:
//   path(distance_matrices) from ch_distance_matrices_for_mds.map{ it [1] }.collect()

//   output:
//   path('*pdf')

//   script:
//   """
//   ${baseDir}/distance_to_MDS_plot.py ${params.four_mN} ${params.n_ind_per_pop} ${distance_matrices}
//   """
// }

// process do_KNN {
//   tag "${variant_set} n${params.n_ind_per_pop}_m${four_mN}_${snp_set}_l${chrom_tag}"
//   publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/${four_mN}/KNN_classification", mode: 'copy'
//   memory '8GB'
// //  time '10m'
//   executor 'local'

  
//   input:
//   tuple variant_set, four_mN, path(distance_matrices) from ch_distance_matrices_for_knn

//   output:
//   path ('KNN*.txt') 
//   path ('') into ch_dummy_plotting_delay // Output goes into a dummy channel that is only used to delay the plotting of all the KNN results across all four_mN runs thus far.

//   script:
//   """
//   ${baseDir}/distance_to_KNN.py ${four_mN} ${params.n_ind_per_pop} ${variant_set} ${params.knn} ${distance_matrices}
//   """
// }

// /* Create a channel that picks up all KNN results from the same chrom_length regardless of four_mN value. 
//   In that channel, mix the dummy delay channel and fiter for unique files to avoid any duplications. */
// ch_for_knn_plotting=ch_dummy_plotting_delay
//           .mix(Channel.fromPath("${baseDir}/../results/n${params.n_ind_per_pop}/${chrom_tag}/*/KNN_classification/KNN_K${params.knn}*.txt"))
//           .unique()
// /*          .dump(tag:"for KNN Plot")*/

// process plot_KNN {
//   tag "n${params.n_ind_per_pop} l${chrom_tag}"
//   publishDir "${baseDir}/../plots/n${params.n_ind_per_pop}/${chrom_tag}/KNN_classification/", mode: 'copy'
//   memory '8GB'
// //  time '10m'
//   executor 'local'


//   input:
//   path(knn_files) from ch_for_knn_plotting.collect()

//   output:
//   path('KNN_summary_plot*pdf')

//   script:
//   """
//   ${baseDir}/plot_KNN.py ${chrom_tag} ${params.n_ind_per_pop} ${knn_files}
//   """
// }

