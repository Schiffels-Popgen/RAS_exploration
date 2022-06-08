#!/usr/bin/env python3

def pop_model_scaled(
    num_replicates=1, 
    sample_sizes=[1,1,1,1,1,1,1,1,1], ## Number of individuals. scaled by times 2 for lineages.
    length=1e6,
    four_mN=1, ## migration per generation is this number times 1.25e-5.
    mutation_rate=1.1e-8,
    recombination_rate=0,
    chrom_name=1
    ):
        Ne=20000
        m = four_mN/(4*Ne)
#         print (m,file=sys.stderr)
        md = 0
#         md = m / sqrt(2) # Diagonal migration, given longer distance in a flat plane.
        ## Matrix of 9 sub-populations
        population_configurations = [
            msprime.PopulationConfiguration(initial_size=Ne, sample_size=2*sample_sizes[i]) for i in range(len(sample_sizes))
        ]
        migration_matrix = [
        ##   P   P   P   P   P   P   P   P   P
        ##   o   o   o   o   o   o   o   o   o
        ##   p   p   p   p   p   p   p   p   p
        ##   1   2   3   4   5   6   7   8   9
            [0,  m,  0,  m,  md, 0,  0,  0,  0],    # Pop1
            [m,  0,  m,  md, m,  md, 0,  0,  0],    # Pop2
            [0,  m,  0,  0,  md, m,  0,  0,  0],    # Pop3
            [m,  md, 0,  0,  m,  0,  m,  md, 0],    # Pop4
            [md, m, md,  m,  0,  m,  md, m,  md],   # Pop5
            [0,  md, m,  0,  m,  0,  0,  md, m],    # Pop6
            [0,  0,  0,  m,  md, 0,  0,  m,  0],    # Pop7
            [0,  0,  0,  md, m,  md, m,  0,  m],    # Pop8
            [0,  0,  0,  0,  md, m,  0,  m,  0]]    # Pop9
        replicates = msprime.simulate(Ne = Ne,
            population_configurations = population_configurations,
            migration_matrix = migration_matrix,
    #         num_replicates = num_replicates,
            mutation_rate = mutation_rate,
            length = length,
            recombination_rate=recombination_rate, random_seed=chrom_name)
        return(replicates)
        

#### MAIN ####

import os
import sys
import msprime
import numpy as np
import pandas as pd
import argparse
from math import ceil as round_up

parser = argparse.ArgumentParser(usage = "%(prog)s [-c chrom_length] [-n chrom_name] [-s sample_size] [-m four_mN] [-o out_dir]" , description = "Simulate data for a structured population whose subpopulations are arranged in a 3x3 grid and have a per-generation migration rate of <four_mN> with their closest non-diagonal neighbours.")
parser._optionals.title = "Available options"
parser.add_argument("-c", type = float, metavar = "chrom_length", required = False, default=1e5, help = "The Length of the chromosome to be simulated.")
parser.add_argument("-n", type = int, metavar = "chrom_name", required = False, default=1, help = "The name of the simulated chromosome. Also functions as the random seed for the simulation.")
parser.add_argument("-s", type = int, metavar = "sample_size", required = False, default=20, help = "The number of samples in each of the four sub-populations.")
parser.add_argument("-m", type = float, metavar = "four_mN", required = False, default=1.0, help = "The scaled per-generation migration rate.")
# parser.add_argument("-o", type = str, metavar = "out_dir", required = False, default="/projects1/MICROSCOPE/rarevar_sim_study/data", help = "The desired output directory within which the freqsum- and eigenstrat-formatted variants will be saved.")
args = parser.parse_args()

# out_base_dir=args.o
four_mN=args.m
sample_sizes=[args.s for _ in range(9)] ## Equal sample sizes for each of the populations.
chrom_name=args.n
chrom_length=args.c

# outdir=out_base_dir+"/{}/{}/".format(np.format_float_scientific(chrom_length, exp_digits=1, trim='-'), four_mN)
outdir='./'

# os.makedirs(outdir, exist_ok=True)


#### SIMULATE DATA ####

trees=pop_model_scaled(
                sample_sizes=sample_sizes,
                four_mN=four_mN,
                mutation_rate=1.1e-8,
                recombination_rate=1e-8,
                length=chrom_length,
                chrom_name=chrom_name
               )

## Extract lineage genotypes from the simulated tree
A=trees.genotype_matrix()
all_positions=np.array([variant.site.position for variant in trees.variants()])

## Pair up lineages from the same population to make diploid individuals

## First create a genotype matrix with half the columns (because individuals should be diploid)
(nrow,ncol) = A.shape
num_inds = np.int(ncol/2)
all_vars = np.zeros((nrow,num_inds), dtype=int)

## Then give each individual a genotype equal to the sum of the allele counts of its lineages.
# print ("IND","H1", "H2", sep="\t")
for ind in range(num_inds):
    haplotype1 = ind*2
    haplotype2 = haplotype1 + 1
#     print (ind, haplotype1, haplotype2, sep="\t")
    all_vars[:, ind] = A [:,haplotype1] + A[:,haplotype2]


#### FILTER DATA ####

## Filter data for commmon, 1M common and rare variants and create separate outputs.
num_lineages=A.shape[1]

## Group definitions:
##  - Common = AF range: [0.4, 0.6]
##  - Rare = AC: >= 10. AF range ![10/num_lineages, -10/num_lineages]

## Convert AF to AC for filtering
## Common Vars
min_ac_common = 0.05 * num_lineages
max_ac_common = 0.95 * num_lineages

### OUTPUT ALL_VARS IN FREQSUM FORMAT INSTEAD. This way we can calculate RAS on rare variants as well as RASCal F3 for all_vars. 
## Rare Vars still calculated for variant_count statistics
max_ac_derived_rare = 10
# min_ac_ancestral_rare = num_lineages - 10 #Rare ancestral alleles (These get filtered out when using RAS anyway)

## Initialise empty filters
common_var_filter = []
rare_var_filter = []

# print(common_vars)
# print (min_ac_rare, max_ac_rare)

## Create filtered arrays with common vs rare genotypes
for s in all_vars.sum(axis=1):
#     print (g , s)
#     print (s)
    if s >= min_ac_common and s <= max_ac_common:
        common_var_filter.append(True)
    else:
        common_var_filter.append(False)
    if s <= max_ac_derived_rare: #or s >= min_ac_ancestral_rare: ## These are filtered by RASCal anyway.
        rare_var_filter.append(True)
    else:
        rare_var_filter.append(False)

common_vars = all_vars[common_var_filter]
common_positions = all_positions[common_var_filter]

rare_vars = all_vars[rare_var_filter]
rare_positions = all_positions[rare_var_filter]

## The first variant count printing will overewrite an existing file, if one exists.
print ('{:<10}'.format(all_vars.shape[0]), "total variants", sep="\t", 
       file=open(outdir+"variant_counts.m{}_chr{}.txt".format(four_mN, chrom_name),'w'))

with open(outdir+"variant_counts.m{}_chr{}.txt".format(four_mN, chrom_name),'a') as f:
    print ('{:<10}'.format(common_vars.shape[0]), "common variants", sep="\t", file=f)
    print ('{:<10}'.format(rare_vars.shape[0]), "rare variants", sep="\t", file=f)
## Twelve-forty variants will need to be generated after all chromosomes of the common variants are concatenated. (should be ok for runtime since the size of these has a maximum value.)


#### EXPORT DATA ####

## Export variants in eigenstrat format
for variant_set_name in ["all_vars", "common_vars"]: ## Twelve_forty will be made later.
    variant_set_genos=abs(eval(variant_set_name)-2) ## Eigenstrat counts Ref (i.e. ancestral) alleles, so the counting needs to be flipped.
    variant_set_positions=eval(variant_set_name.replace("vars","positions"))
    ## Export simulated eigenstrat geno file for testing with qp3Pop
    np.savetxt(outdir+variant_set_name+"_m{}_chr{}.geno".format(four_mN, chrom_name), variant_set_genos, delimiter='', fmt="%d")

    ## Export .snp file to go with it.
    #### Resulting snp file is unsorted, so it needs sorting within bash to function.
    with open(outdir+variant_set_name+"_m{}_chr{}.snp".format(four_mN, chrom_name), 'w') as f:
        for pos in variant_set_positions:
            print("{}_{}".format(chrom_name, round_up(pos)), str(chrom_name), "{:.5f}".format(pos/chrom_length), round_up(pos), "A", "G", sep="\t", file=f)

    ## Export .ind file
    with open(outdir+variant_set_name+"_m{}_chr{}.ind".format(four_mN, chrom_name), "w") as f: 
        ## pop defined as msprime populations from 0 to 2x#inds (so the entire list since individuals are diploid), 
        ## but considering every second population in the list, since two lineages are merged to form an individual.
        for ind,pop in enumerate(trees.tables.nodes.population[0:sum(sample_sizes)*2:2]):
            print("ind"+str(ind), "U", "Pop"+str(pop), sep="\t", file=f)

## Create a freqsum file out of the rare variants
#### Ind "Ref" is created with all genotypes being ancestral ("0")
freqsum_output=outdir+"/all_vars_m{}_chr{}.freqsum".format(four_mN, chrom_name)
# freqsum_output=outdir+"/rare_vars_m{}_chr{}.freqsum".format(four_mN, chrom_name)
with open(freqsum_output, 'w') as f:
    print ("#CHROM","POS","REF","ALT",*["ind"+str(x)+"(2)" for x in range(num_inds)],"Ref(2)", sep="\t", file=f)
    for pos,genos in zip(all_positions, all_vars):
        print(str(chrom_name), round_up(pos), "A", "G", *genos, "0", sep="\t", file=f)
## Older generation of rarevars only freqsum. Now outputting all_vars freqsum instead.
    # for pos,genos in zip(rare_positions, rare_vars):
    #     print(str(chrom_name), round_up(pos), "A", "G", *genos, "0", sep="\t", file=f)