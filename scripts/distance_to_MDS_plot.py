#!/usr/bin/env python3

## A function to carry out PCoA and plot the results
## Function to read in matrices. Returns a dictionary of the form { "snp_set" : distance_matrix}
def read_matrices(inputs): ## inputs should be a list of file names
    dm={}
    for file in inputs:
        snp_set=re.sub(r'_[a-zA-Z_0-9.]*$','', file) ## Infer snp_set from file name
        input_matrix = pd.read_csv(file, sep="\t", header=None)
        # np.fill_diagonal(input_matrix.values, 0) ## Matrices should already have 0 diagonals.
        dm[snp_set]=input_matrix
    return(dm)


## A function to carry out PCoA and plot the results
def do_mds(dm, m, population_colours, num_dim = 2, num_inds_per_pop = 20):
    ## dm is a dictionary where the keys are the snp sets and values are the distance matrix pd dataframes. 
    x_2d={}
    output_plot="MDS_m{}.pdf".format(m)
    f, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4)

    f.set_figheight(7.5)
    f.set_figwidth(4*7.5)
    plot_column={"all":ax1, "common":ax2, "1240k":ax3, "rare":ax4}

    for snp_set,matrix in dm.items():
        x_2d[snp_set]=pcoa(matrix, number_of_dimensions = num_dim)
        add_subplot(plot_column[snp_set], x_2d[snp_set], snp_set, m, population_colours)
    legend_elements = [Line2D([0], [0], marker='o', color='w', label="Pop"+str(pop),
                              markerfacecolor=population_colours["Pop"+str(pop)], markersize=10) for pop in range(9)]
    lgd=f.legend(handles=legend_elements, loc='right',bbox_to_anchor=(1.00, 0.5), prop={"size":14})
    f.savefig(output_plot, bbox_extra_artists=(lgd,), bbox_inches='tight')

## Function to add the given data to a specific subplot with the correct annotations
def add_subplot(ax, x_2d, snp_set, m, population_colours):
#     ax.set(aspect='equal')
    ax.axis('equal')
    ax.set_xlabel("Coordinate 1")
    ax.set_ylabel("Coordinate 2")
    ax.scatter(x_2d.samples["PC1"],x_2d.samples["PC2"], 
                c=np.repeat(list(population_colours.values()), 10),
               alpha=0.5)
    ax.title.set_text("{} Variants (m={})".format(snp_set.capitalize(), m))
    

##### MAIN #####

import sys
from skbio.stats.ordination import pcoa ## Classical MDS
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import re
import numpy as np
## Import population colours from auxiliary file to make changing colours easy in the future.
from colour_definitions import population_colours

args=sys.argv

if "-h" in args or "--help" in args or len(sys.argv) == 1:
  print("""
  A small script to carry out PCoA on a set of input distance matrices and plot coordinate 1 vs 2 in a 4x1 compound plot for each variant. 
  All arguments are positional. Multiple distance matrices can be provided one after another.
  Only one distance matrix is expected for the rare variants.
  
  Usage:
  distance_to_MDS_plot.py m n_ind_per_pop input_matrix [input_matrix2 ..]
  
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

# print("m={}\tn={}".format(m,n_ind_per_pop))
# print("inputs={}".format(input_matrices))
# sys.exit(0)

dm = read_matrices(input_matrices)
do_mds(dm, m, population_colours, num_dim = 2, num_inds_per_pop = n_ind_per_pop)