#!/usr/bin/env bash

if [[ ${1} == "-h" || ${1} == "--help" || ${1} == '' ]]; then
  echo "This script will add a Reference individual (with all 2 genotypes) to the simulated genotype datasets created by simulation.py,"
  echo "  create an ungrouped ind file for that dataset (.indEach), and sort the genotype and snp files by variant chromosome position." 
  echo ''
  echo "Usage:"
  echo "    $0 chrom_name input_dir"
  exit 0
fi

chrom_name=$1
input_dir=$2
four_mN=$(echo ${2%/} | rev | cut -f 1 -d "/" | rev) ## Infer four_mN from input dir


for var_set in "all_vars" "common_vars"; do

    ## Add fake individual Ref to .ind file, and create indEach file with each individual belonging to own population
    echo -e "Ref\tU\tRef" >>${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.ind
    awk -v OFS="\t" '$3=$1' ${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.ind >${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.indEach

    ## Infer number of SNPs
    num_snps=$(wc -l ${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.snp| cut -d " " -f1)

    ## Create genotypes for "Ref", consisting of ancestral allele across all SNPs ("2")
    for i in `seq 1 1 ${num_snps}`; do echo "2"; done | paste -d '' ${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.geno - >${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.geno2
    mv ${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.geno2 ${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.geno

    ## Sort the .snp file. Keep genotype order in geno as well. They tend to already be sorted, but sometimes that doesnt seem to be the case, so best to be safe.
    paste -d " "   ${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.snp ${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.geno | sort -nk4 >${input_dir}/temp_${var_set}_m${four_mN}_chr${chrom_name}.snpgeno
    cut -d " " -f1 ${input_dir}/temp_${var_set}_m${four_mN}_chr${chrom_name}.snpgeno >${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.snp2
    cut -d " " -f2 ${input_dir}/temp_${var_set}_m${four_mN}_chr${chrom_name}.snpgeno >${input_dir}/${var_set}_m${four_mN}_chr${chrom_name}.geno2
    ## Remove temp snpgeno file.
    rm ${input_dir}/temp_${var_set}_m${four_mN}_chr${chrom_name}.snpgeno
    
    ### For small chromosome lengths (i.e. during basic testing) sometimes two snps will get rounded to the same position, and sorting will swap around the order of genotypes randomly. Remains to be seen if this still occurrs when the lengths are larger. This happens because msprime gives unique float positions to sites, and these need to be rounded up to make them bp positions.
done
