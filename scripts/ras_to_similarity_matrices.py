#!/usr/bin/env python3
##  A function to add the RAS and NSites of a ras log to a matrix of ras results.
##  If the ras result matrix is empty, returns the given ras log.
def sum_ras_logs(summed_logs, ras_log, max_af):
    ## Read in data while ignoring the first 6 lines in the file (header)
    data = pd.read_csv(ras_log, sep="\t", skiprows=6, header=None, names=["IndA","IndB","RAS","NSites","ras","rasErr","AC"])
    ## Delete rows where the measurement is Outgroup F3 (i.e. keep only ras)
    data = data[data["AC"]!="Outgroup F3"]
    ## Delete the total RAS rows, since these need to be recalculated anyway.
    data = data[data["AC"]!="Total [2,{}]".format(max_af)]
    ##Now the AC column can be converted to int
    data["AC"] = data["AC"].astype(int)
    
    ## To save on memory, drop the unused ras and rasErr columns.
    ## Ras will be recalculated from the RAS and NSites columns after all chromosomes are compiled. 
    data = data.drop(["ras","rasErr"],axis=1)
    ## If summed_logs is empty, return the parsed ras_log contents
    if summed_logs.empty:
      return(data)
    ## Otherwise, add the "RAS" and "NSites" of the two matrices and return the new matrix.
    else:
      summed_logs["RAS"] = summed_logs["RAS"] + data["RAS"]
      summed_logs["NSites"] = summed_logs["NSites"] + data["NSites"]
      return(summed_logs)

## Function to convert the summed ras logs into a ras similarity matrix.
def convert_to_similarity_matrix(summed_logs, max_af):
    ## First, compute ras from the RAS and NSites columns
    summed_logs["ras"] = summed_logs["RAS"] / summed_logs["NSites"]
    ## These columns are no longer needed
    summed_logs = summed_logs.drop(["RAS","NSites"], axis=1)
    
    summed_ras_by_ac={}
    for m in range(2, max_af+1):
      ## Then group by IndA/IndB combination and sum ras up to each AC
      ## Keep only lines up to desired AC, group by IndA/IndB combination and sum those lines to get summed ras (and AC which is nonsense).
      ## rename_axis and reset_index reformat the resulting dataframe to include IndA and IndB as actual columns and not rownames.
      grouped_filtered = summed_logs[ summed_logs["AC"] <= m ].groupby(["IndA","IndB"]).sum().rename_axis(["IndA","IndB"]).reset_index()
      
      ## Infer number of individuals
      num_inds = len(grouped_filtered["IndA"].unique())
      
      ## Pivot grouped_filtered into a similarity matrix
      grouped_filtered = grouped_filtered.pivot(index="IndA", columns="IndB", values="ras")
      ## The matrix then needs to be reindexed so it is in the desired order (by individual number instead of alphanumerically)
      desired_order = ["ind"+str(i) for i in range(num_inds)]
      grouped_filtered = grouped_filtered.reindex(index=desired_order, columns=desired_order)
      
      ## Set diagonal (within individual) ras to 1?
      # np.fill_diagonal(filtered_data.values, 1)
      
      summed_ras_by_ac[m]=grouped_filtered
    return(summed_ras_by_ac)

def export_ras_matrices(summed_ras_by_ac):
    for af,matrix in summed_ras_by_ac.items():
      matrix.to_csv("rare_similarity_matrix.ac{}.txt".format(af), sep='\t', header=False, index=False, float_format='%.9f')

import pandas as pd
import argparse
import sys

args=sys.argv[1:]

if "-h" in args or "--help" in args or len(sys.argv) == 1:
  print("""
  A small script to compile cumulative AC similarity matrices out of a provided set of ras output files.
  All arguments are positional. Multiple log files can be provided one after another.
  
  Usage:
  ras_to_similarity_matrices.py <Max_ras_ac> <ras output files>
  
  Options:
    -h,--help   Print this text and exit.
  """)
  sys.exit(0)


try:
  max_af=int(args[0])
except ValueError:
  print("""
  ERROR: The first argument should be an integer equal to the maximum allele count used for RASCal (-M).
  Execution halted.
  """)

ras_logs=args[1:]

## Sum RAS and NSites across all chromosomes
summed_logs=pd.DataFrame()
for log in ras_logs:
  summed_logs = sum_ras_logs(summed_logs, log, max_af)

## Calcuate ras from summed RAS and NSites, then convert to a similarity matrix
summed_ras_by_ac = convert_to_similarity_matrix(summed_logs, max_af)

## Export matrices into tsv format
export_ras_matrices(summed_ras_by_ac)