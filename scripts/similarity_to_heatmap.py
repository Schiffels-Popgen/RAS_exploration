#!/usr/bin/env python3

## Function to read in similarity matrices. Returns a dictionary of the form { "snp_set" : distance_matrix}
def read_matrices(inputs): ## inputs should be a list of file names. Returns a dictionary 
    sm={}
    for file in inputs:
        snp_set=re.sub(r'_[a-zA-Z_0-9.]*$','', file) ## Infer snp_set from file name
        input_matrix = pd.read_csv(file, sep="\t", header=None)
        np.fill_diagonal(input_matrix.values, None) ## Add NA in the diagonal. Otherwise all other colors are not bright enough.
        sm[snp_set]=input_matrix
    return(sm)

## Function to plot a heatmap from a set of similarity matrices.
def plot_heatmaps(sm, m, n_ind_per_pop, population_colours):
    ## sm is a dictionary where the keys are the snp sets and values are the distance matrix pd dataframes. 
    n_ind_per_pop = n_ind_per_pop
    num_inds = 9 * n_ind_per_pop
    f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(ncols=5)
    
    f.set_figheight(10)
    f.set_figwidth(5*10)
    plot_column={"all":ax1, "all.rascal":ax2, "common":ax3, "1240k":ax4, "rare":ax5}
    plt.tight_layout(rect=(0,0,1,0.975))
    
    for snp_set,matrix in sm.items():
        plt.rcParams.update({'font.size': 20})
        add_subplot(plot_column[snp_set], sm[snp_set], snp_set, m, population_colours, num_inds, n_ind_per_pop)
    f.savefig('Heatmap.m{}.pdf'.format(m))

## Helper function to add a subplot to the image
def add_subplot(ax, data, snp_set, m, population_colours, num_inds, n_ind_per_pop):
    ax.set_aspect('equal')
    ax.imshow(data)
    for i in ax.spines:
        ax.spines[i].set_visible(False)
    ax.set_xticks(np.arange(0,num_inds, n_ind_per_pop)-0.5, minor=True)
    ax.set_yticks(np.arange(0,num_inds, n_ind_per_pop)-0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.xaxis.set_tick_params(width=1.5, length=6)
    ax.yaxis.set_tick_params(width=1.5, length=5)

    ax.title.set_text("{} Variants (m={})".format(snp_set.capitalize(), m))

##### MAIN #####

import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd
from colour_definitions import population_colours

args=sys.argv

if "-h" in args or "--help" in args or len(sys.argv) == 1:
  print("""
  A small script to plot a set of input distance matrices as a heatmap, and generate a 4x1 compound plot for each variant. 
  All arguments are positional. Multiple distance matrices can be provided one after another.
  Only one similarity matrix is expected for the rare variants.
  
  Usage:
  similarity_to_heatmap.py m n_ind_per_pop input_matrix [input_matrix2 ..]
  
  Options:
    -h,--help   Print this text and exit.
  """)
  sys.exit(0)

try:
  m=float(args[1])
except ValueError:
  print("""
  ERROR: The first argument should be a float equal to the scaled migration rate used for the simulation.
  Execution halted.
  """)

try:
  n_ind_per_pop=int(args[2])
except ValueError:
  print("""
  ERROR: The second argument should be an integer equal to the number of sampled individuals in each of the 9 simulated populations.
  Execution halted.
  """)
input_matrices=args[3:]

sm=read_matrices(input_matrices)
plot_heatmaps(sm, m, n_ind_per_pop, population_colours)