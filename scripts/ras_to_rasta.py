#!/usr/bin/env python3

## Function to infer chromosome number from input file name
def infer_chrom_number(input_file):
  ## Infer chromosome number form file name
  chrom = input_file.split("_")[-1].split(".")[0].lstrip("chr")
  return (chrom)

## Function that reads an input ras output and converts the ind-ind sharing matrix to ind-pop ras.
##     Returns a dictionary with structure {AC : ind-pop_ras_matrix}
def read_in_ras(input_file, max_af):
    ## Read in data
    ras_data = pd.read_csv(input_file, sep="\t", skiprows=6, header=None, 
                           names=["IndA","IndB","RAS","NSites","ras","rasErr","AC"])
    ## Delete rows where the measurement is Outgroup F3 (i.e. keep only ras)
    ras_data = ras_data[ras_data["AC"]!="Outgroup F3"]
    ## Delete the total RAS rows, since these need to be recalculated anyway.
    ras_data = ras_data[ras_data["AC"]!="Total [2,{}]".format(max_af)]
    ##Now the AC column can be converted to int
    ras_data["AC"] = ras_data["AC"].astype(int)
    ## To save on memory, drop the unused ras and rasErr columns.
    ## Ras will be recalculated from the RAS and NSites columns after all chromosomes are compiled. 
    ras_data = ras_data.drop(["ras","rasErr"],axis=1)

    if len(ras_data.NSites.unique()) != 1:
        raise ValueError("Non-equal number of RAS sites within each chromosome.")
    else:
        n_sites=ras_data.NSites.unique()[0]

    ## Infer number of individuals
    ind_list = ras_data["IndA"].unique()
    num_inds = len(ind_list)
    inds_per_pop = int(num_inds/9)
    ## Use number of individuals to get a list of inds per population
    ## This will be used to check if an individual is within a population to adjust population sizes in RAS calculation
    pops = {}
    for p in range(9):
        pops["Pop"+str(p)] = ind_list[p*inds_per_pop:p*inds_per_pop+inds_per_pop]

    ## Remove NSites column. no longer needed
    ras_data = ras_data.drop(["NSites"],axis=1)

    ## Split ras and sum up to each AC. 
    ras_data_by_ac = {}
    for m in range(2, max_af+1):
        sum_ras = ras_data [ras_data["AC"] <= m].groupby(["IndA","IndB"]).sum().rename_axis(["IndA","IndB"]).reset_index()

        ## Convert data to wide format, set within individual sharing to 0. 
        sum_ras = pd.pivot_table(sum_ras, columns="IndB", values="RAS",  index="IndA")
        sum_ras = sum_ras[["ind"+str(i) for i in range(90)]].reindex(["ind"+str(i) for i in range(90)])
        np.fill_diagonal(sum_ras.values,0)
        ## Sum individual columns into populations based on individual list for each population
        for p,inds in pops.items():
            sum_ras[str(p)] = sum_ras[inds].sum(axis=1)
        sum_ras = sum_ras.iloc[:,-9:] ## Keep only Population columns
        ## apply edits the dataframe in place
        sum_ras.apply(normalise_pop_sizes, axis=1, result_type="expand", pops=pops)

        ras_data_by_ac[m]=sum_ras
    return(ras_data_by_ac,pops,inds_per_pop, n_sites)

## Function to calculate group sharing profiles for each group
def calculate_group_sharing(ras_data_by_ac,pops, inds_per_pop):
    grouped_ras_by_ac={}
    for m in range(2, max_af+1):
        sum_ras = ras_data_by_ac[m]
        ## Add a column with the population assignment of each individual
        grouped_ras = sum_ras.apply(assign_pop_label, axis=1, pops=pops)
        ## Sum columns by population label and normalise by number of individuals per population
        grouped_ras = grouped_ras.groupby("Group").sum()/inds_per_pop
        grouped_ras_by_ac[m]=grouped_ras
    return (grouped_ras_by_ac)

## Function to convert RAS table to RASTA by calculating group and individual differences.
def convert_ras_to_rasta(ras_data_by_ac, grouped_ras_by_ac, pops, inds_per_pop, n_sites, blockSizes_dict):
    rasta_by_ac = {}
    for m,ras_data in ras_data_by_ac.items():
        rasta_table = ras_data.apply(
            assign_pop_label,
            axis=1, 
            pops=pops
            ).apply(
            calculate_rasta, 
            axis=1, 
            grouped_ras=grouped_ras_by_ac[m], 
            inds_per_pop=inds_per_pop
           )
        rasta_table = pivot_rasta_table(rasta_table, n_sites, chrom, blockSizes_dict)
        rasta_by_ac[m] = rasta_table
    return (rasta_by_ac)

## Function that takes the individual and Group RAS profiles and returns RASTA (indX, truePop; otherPops, Ref)
## Should be applied rowwise to the labelled ras_data 
def calculate_rasta(labelled_ras_row, grouped_ras, inds_per_pop):
    true_pop = labelled_ras_row["Group"]
    group_sharing = grouped_ras.loc[true_pop]
    ## In order to get the sharing that the group would have without that 1 individual, we need to multiply the sharing by inds_per_pop, subtract the individual and multiply by n-1
    row_rasta = (labelled_ras_row - (inds_per_pop * group_sharing))/(inds_per_pop-1)
    return(row_rasta.drop("Group"))

## Make RASTA table into long format
def pivot_rasta_table(rasta_df, n_sites, chrom, blockSizes_dict):
    ## Add pop groups and make IndA into a column
    rasta_df = rasta_df.apply(assign_pop_label, axis=1, pops=pops).reset_index()
    ## Convert to long format data
    rasta_df = pd.melt(rasta_df, id_vars=["IndA","Group"], value_name="RASTA",
                       value_vars=["Pop0", "Pop1","Pop2","Pop3","Pop4","Pop5","Pop6","Pop7","Pop8"])
    ## Add uniform nSites & Ref column, rename columns so f4 equivalence is clearer
    rasta_df["nSites"] = n_sites
    rasta_df["Chrom"] = chrom
    rasta_df["blockSize"] = blockSizes_dict[chrom]
    rasta_df["Ref"] = "Ref"
    rasta_df = rasta_df.rename(columns={"IndA":"Ind","Group":"PopA","IndB":"PopB"})
    ## Return the dataframe with column order tweaked to make more legible
    return(rasta_df[["Ind","PopA","PopB","Ref","RASTA","nSites", "Chrom", "blockSize"]])

## Function to normalise the sharing with each population by the correct size of the populations without the target individual
def normalise_pop_sizes(row, pops):
    for pop in pops.keys():
        row[pop]= row[pop]/len(np.setdiff1d(pops[pop],row.name))
    return (row)

## Function to add a "Group" column with the true population of that individual.
def assign_pop_label(row, pops):
    for pop,inds in pops.items():
        if row.name in inds:
            row["Group"] = pop
    return row

## Jackknife function from RASUtils
  ## blockValues = list of RASTA value per chromosome
  ## totalObservations = list of n_sites per chromosome
  ## blockSizes = list of Freqsum lengths per chromosome
def getJackknife(blockValues, totalObservations, blockSizes, skipJackknife):
    thetaminus=[0 for x in range(len(blockSizes))]
    sum1=0
    sum2=0
    jackknifeStdErr=0
    if sum(totalObservations)==0:
        ## If no observations were made, the rare allele sharing rate is 0.
        thetahat=0
    else:
        ## thetahat is normalised by number of observations
        thetahat=sum(blockValues)/sum(totalObservations)

    for c in range(len(blockSizes)):
        if totalObservations[c] == sum(totalObservations):
            ## If all rare alleles are on a single chunk, the ras/site without that chunk is 0.
            thetaminus[c]=0
        else:
            ## thetaminus is normalised by number of observations
            thetaminus[c]=( (sum(blockValues)-blockValues[c]) / (sum(totalObservations)-totalObservations[c]) )
        
        ## Jackknife estimator calculation
        sum1+=thetahat-thetaminus[c]
        sum2+=(blockSizes[c]*thetaminus[c])/sum(blockSizes)
    jackknifeEstimator=sum1+sum2
    
    ## Standard error calculation
    if not skipJackknife:
        for c in range(len(blockSizes)):
            hj=sum(blockSizes)/blockSizes[c]
            pseudoval = (hj*thetahat)-((hj-1)*thetaminus[c])
            jackknifeStdErr+=(1/len(blockSizes))*(((pseudoval-jackknifeEstimator)**2)/(hj-1))
    return (jackknifeEstimator,jackknifeStdErr)

def calculate_ZScore(thetaJ, sigma2):
  try:
      return (thetaJ/sigma2)
  except ZeroDivisionError:
      return (0)
  

##### MAIN #####

import sys
import numpy as np
import pandas as pd

args=sys.argv

if "-h" in args or "--help" in args or len(sys.argv) == 1:
  print("""
  A small script to calculate RASTA from RAS.
  All arguments are positional. Multiple distance matrices can be provided one after another.
  Output file names are hard-coded.
  
  Usage:
  ras_to_rasta.py max_af blockSizes input_ras [input_ras2 ..]
  
  Options:
    -h,--help   Print this text and exit.
  """)
  sys.exit(0)

max_af = int(args[1])
blockSizes_file = args[2]
input_files = args[3:]

## Read in blockSizes from file
blockSizes_dict = {}
for line in open(blockSizes_file, 'r'):
  fields = line.strip().split()
  chrom = fields[0].split("_")[-1].split(".")[0].lstrip("chr")
  blockSizes_dict[chrom] = int(fields[1])-1

## Initialise final dict
rasta_results_per_chrom={}

for input_file in input_files:
  ## Infer Chrom number from file name
  chrom = infer_chrom_number(input_file)
  
  ## Read in RAS and infer pop definitions, number of inds per pop and n_sites for this chrom
  (results, pops, inds_per_pop, n_sites) = read_in_ras(input_file, max_af)
  grouped_results = calculate_group_sharing(results, pops, inds_per_pop)
  
  ## Save converted RASTA to dictionary of form {chrom : { AC : RASTA } }
  rasta_results_per_chrom[chrom] = convert_ras_to_rasta(results, grouped_results, pops, inds_per_pop, n_sites, blockSizes_dict)

## Concatenate all RASTAs together per AF
concatenated_rasta_by_ac = {}
dfs_to_concat = []

## rasta_table is { ac : RASTA } for each chromosome
for ac in range (2,max_af+1):
  for rasta_table in rasta_results_per_chrom.values():
    ## Make a list with all the dfs that should be concatenated. More efficient than appending to a df in a for loop.
    dfs_to_concat.append(rasta_table[ac])
  ## Concatenate all rasta dfs across all chroms for each ac
  concatenated_rasta_by_ac[ac] = pd.concat(dfs_to_concat)

## Apply Jackknife and print results
for ac in range(2,max_af+1):
  out_file = open("rasta.m{}.txt".format(ac), 'w')
  ## Print output header
  print ("Ind", "PopA", "PopB", "Ref", "RASTA", "stderr", "ZScore", "nSites", sep="\t", file=out_file)
  for idx in range(inds_per_pop * 9):
    ind = "ind"+str(idx)
    population_labels= ["Pop"+str(_) for _ in range(9)]
    for popA in population_labels:
      for popB in population_labels:
        subset = concatenated_rasta_by_ac[ac].query(
          '(Ind == "{}") & \
          (PopA == "{}") & \
          (PopB == "{}")'.format(ind, popA, popB)
        )
        # print(subset)
        # print(list(subset.RASTA), list(subset.nSites), list(subset.blockSize))
        thetaJ,sigma2 = getJackknife(list(subset.RASTA), list(subset.nSites), list(subset.blockSize), False)
        ZScore = calculate_ZScore(thetaJ,sigma2)
        print (ind, popA, popB, "Ref", thetaJ, sigma2, ZScore, sum(subset.nSites), sep="\t", file=out_file)