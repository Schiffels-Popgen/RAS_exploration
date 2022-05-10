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

// Parameters to refer to within multiple processes. 
num_inds = params.n_ind_per_pop * 9 - 1
package_dir = "${baseDir}/../data/n${params.n_ind_per_pop}/${params.chrom_length}/${params.four_mN}/poseidon/"

// Find simulated data and create input channels from params.
  // This should snsure execution only when input files change.
ch_input_datasets=Channel.fromList(["all_vars", "common_vars", "twelve_forty_vars", "rare_vars"])
  .map{
    it ->
      def x = it == "rare_vars" ? "all_vars" : it //Pick up all_vars dataset for rare variants.
      def bed = file("${package_dir}/"+x+"/"+x+".bed")
      def bim = file("${package_dir}/"+x+"/"+x+".bim")
      def fam = file("${package_dir}/"+x+"/"+x+".fam")

    [it, bed, bim, fam]
  }
  .branch{
    rare: it[0] == "rare_vars"
    common: true
  }
  .into{ 
    ch_input_for_xerxes_ras; ch_input_xerxes_pairwise_ras; ch_input_for_xerxes_f3; ch_input_for_xerxes_f4
  }

process xerxes_pairwise_ras {
  tag "n${params.n_ind_per_pop}_m${params.four_mN}_l${params.chrom_length}"
  publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${params.chrom_length}/${params.four_mN}/xerxes_ras", mode: 'copy'
  memory '8GB'
  cpus 1

  input:
  tuple variant_set, file(bed), file(bim), file(fam) from ch_input_xerxes_pairwise_ras.rare

  output:
  tuple variant_set, path("pairwise_*.out") into (ch_xerxes_pairwise_ras_output_for_matrix)
  file "pairwise_popConfigFile.txt"
  file "pairwise_blockTableFile.txt"

  script:

  """
  ## Create population definitions file. Each pop consists of all but the last individual. The first individual of a pop is used as a left pop.
  echo "groupDefs:" >pairwise_popConfigFile.txt
  leftpops=()
  rightpops=()
  n_inds=$((${params.n_ind_per_pop} * 9 - 1))

  ## Define rights and lefts
  echo "popLefts:" >>pairwise_popConfigFile.txt
  for idx in \$(seq 0 1 \${n_inds}); do
    echo "  - <ind\${idx}>" >>pairwise_popConfigFile.txt
  done

  echo "popRights:" >>pairwise_popConfigFile.txt
  for idx in \$(seq 0 1 \${n_inds}); do
    echo "  - <ind\${idx}>" >>pairwise_popConfigFile.txt
  done

  echo "outgroup: <Ref>" >>pairwise_popConfigFile.txt

  ## Run ras
  ${params.poseidon_exec_dir}/xerxes ras -d ${package_dir}/all_vars \
    --popConfigFile pairwise_popConfigFile.txt \
    --minAC 2 \
    --maxAC ${params.max_ras_ac} \
    -j CHR \
    -f pairwise_ras_table.out \
    --blockTableFile pairwise_blockTableFile.txt
  """
}

process xerxes_ras {
  tag "n${params.n_ind_per_pop}_m${params.four_mN}_l${params.chrom_length}"
  publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${params.chrom_length}/${params.four_mN}/xerxes_ras", mode: 'copy'
  memory '8GB'
  cpus 1

  input:
  tuple variant_set, file(bed), file(bim), file(fam) from ch_input_for_xerxes_ras.rare

  output:
  tuple variant_set, path("*.out") into (ch_xerxes_ras_output_for_matrix, ch_xerxes_ras_for_rasta)
  file "popConfigFile.txt"
  file "blockTableFile.txt"

  script: 
  """
  ## Create population definitions file. Each pop consists of all but the last individual. The first individual of a pop is used as a left pop.
  echo "groupDefs:" >popConfigFile.txt
  leftpops=()
  rightpops=()
  for pop in {0..8}; do
    excluded_ind='-<ind'\$((${params.n_ind_per_pop} * \${pop}))'>'
    echo -e "  Pop\${pop}Rest: Pop\${pop},\${excluded_ind}" >>popConfigFile.txt
    leftpops+=(\${excluded_ind#"-"}) ## Add excluded ind to lefts without the leading "-"
    rightpops+=(Pop\${pop}Rest) ## Add created population into rights
  done

  ## Define rights and lefts
  echo "popLefts:" >>popConfigFile.txt
  for left in \${leftpops[@]}; do
    echo "  - \${left}" >>popConfigFile.txt
  done
  for right in \${rightpops[@]}; do
    echo "  - \${right}" >>popConfigFile.txt
  done

  echo "popRights:" >>popConfigFile.txt
  for right in \${rightpops[@]}; do
    echo "  - \${right}" >>popConfigFile.txt
  done

  echo "outgroup: <Ref>" >>popConfigFile.txt

  ## Run ras
  ${params.poseidon_exec_dir}/xerxes ras -d ${package_dir}/all_vars \
    --popConfigFile popConfigFile.txt \
    --minAC 2 \
    --maxAC ${params.max_ras_ac} \
    -j CHR \
    -f ras_table.out \
    --blockTableFile blockTableFile.txt
  """
}

process xerxes_f3 {
  tag "${variant_set} n${params.n_ind_per_pop}_m${params.four_mN}_l${params.chrom_length}"
  publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${params.chrom_length}/${params.four_mN}/xerxes_f3", mode: 'copy'
  memory '8GB'
  cpus 1

  input:
  tuple variant_set, file(bed), file(bim), file(fam) from ch_input_for_xerxes_f3.common

  output:
  tuple variant_set, path("*.out") into (ch_xerxes_fstats_output_for_matrix, ch_xerxes_fstats_for_rasta)

  script:
  """
  ## Create statFile
  ## Create overall poplist
  for i in \$(seq 0 1 ${num_inds}); do
    lower_bound=\$(expr \${i} + 1)
    for j in \$(seq \${lower_bound} 1 ${num_inds}); do
      echo -e "F3vanilla(<ind\${i}>, <ind\${j}>, Ref)"
    done
  done >pairwise_poplist.txt

  ## Remove the last byte (i.e. the final newline character) from the statFile
  truncate -s -1 pairwise_poplist.txt

  ## Run f3
  ${params.poseidon_exec_dir}/xerxes fstats -d ${package_dir}/${variant_set} \
    --statFile pairwise_poplist.txt \
    -j CHR \
    -f f3_table_${variant_set}.out
  """
}

process xerxes_f4 {
  tag "${variant_set} n${params.n_ind_per_pop}_m${params.four_mN}_l${params.chrom_length}"
  publishDir "${baseDir}/../results/n${params.n_ind_per_pop}/${params.chrom_length}/${params.four_mN}/xerxes_f3", mode: 'copy'
  memory '8GB'
  cpus 1

  input:
  tuple variant_set, file(bed), file(bim), file(fam) from ch_input_for_xerxes_f4.common

  output:
  tuple variant_set, path("*.out") into (ch_xerxes_f4_output)

  script:
  """
  ## Create statFile
  for i in \$(seq 0 1 ${num_inds}); do
    for b in \$(seq 0 1 9); do
      for c in \$(seq 0 1 9); do
        if [[ ! "\${b}" == "\${c}" ]]; then
          echo -e "F4(<ind\${i}>, Pop\${b}, Pop\${c}, Ref)"
        fi
      done
    done
  done >pairwise_poplist.txt

  ## Remove the last byte (i.e. the final newline character) from the statFile
  truncate -s -1 pairwise_poplist.txt

  ## Run f4
  ${params.poseidon_exec_dir}/xerxes fstats -d ${package_dir}/${variant_set} \
    --statFile pairwise_poplist.txt \
    -j CHR \
    -f f4_table_${variant_set}.out
  """
}