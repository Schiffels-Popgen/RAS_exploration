#!/usr/bin/env python3 

## Function to parse the input data from multiple KNN logfiles into a single compound table.
def parse_knn_output(input_knn_file, summary_table):
    data=pd.read_csv(input_knn_file, sep="\t")
    ## Extract four_mN and snp_set from the input file (and flatten to single value).
    four_mN=data["four_mN"].unique()[0]
    snp_set=data["variant_set"].unique()[0]
    ## Create disctionary with information. pd.concat will then use this dictionary to add entries into growing dataframe in the respective columns.
    df={'snp_set':snp_set,
        "m" : four_mN,
        "success_rate":[sum(data["correct_guess_proportion"])/data.shape[0]]
       }
    df = pd.concat([pd.DataFrame(df),summary_table], ignore_index=True)
    return (df)

## Function to create a compound KMM plot from the parsed input data.
def plot_KNN (summary_table, snp_set_colours, chrom_length):
    ax = summary_table.pivot(index="m", columns="snp_set", values="success_rate").plot.line(
      color=snp_set_colours, figsize=(15,10))
    ax.legend(title="Variant set")
    ax.set_xscale('log')
    ax.set_ylim(-0.05,1.05)
    fig = ax.get_figure()
    fig.savefig("KNN_summary_plot.l{}.pdf".format(chrom_length))

##### MAIN #####
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from colour_definitions import snp_set_colours

args=sys.argv

if "-h" in args or "--help" in args or len(sys.argv) == 1:
  print("""
  A script to plot the results of multiple KNN runs of the same chromosome length.
  All arguments are positional. Multiple KNN logfiles can be provided one after another.
  
  Usage:
  plot_KNN.py chrom_length input_matrix [input_matrix2 ..]
  
  Options:
    -h,--help   Print this text and exit.
  """)
  sys.exit(0)

try:
  chrom_length=float(args[1])
except ValueError:
  print("""
  ERROR: The first argument should be an integer equal to the length of the simulated choromosomes in bp.
  Execution halted.
  """)

input_files=sys.argv[2:]

summary_table = pd.DataFrame(columns=['snp_set', 'm', 'success_rate'])

for file in input_files:
    summary_table=parse_knn_output(file, summary_table)

plot_KNN (summary_table, snp_set_colours, chrom_length)
