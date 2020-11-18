#!/usr/bin/env python3

import numpy as np
import argparse
import sys

args=sys.argv

if "-h" in args or "--help" in args or len(sys.argv) == 1:
  print("""
  A small script to convert input similarity matrices into distance matrices by returning the 1-matrix.
  All arguments are positional. Multiple log files can be provided one after another.
  
  Usage:
  similarity_to_distance.py <variant_set> <similarity matrices>
  
  Options:
    -h,--help   Print this text and exit.
  """)
  sys.exit(0)


variant_set=args[1]
input_files=args[2:]

## To get a dissimilarity matrix, simply return the 1-X for X in a similarity matrix 
for f in input_files:
  ## Output file name. If the run is on rare variants, infer the max ac of the ras matrix from the file name.
  if variant_set == "rare":
    max_ac=f.split('.')[-2][-1] ## AC is the last character before the last '.'
    output_file_name='{}_distance_matrix.ac{}.txt'.format(variant_set, max_ac)
  else:
    output_file_name='{}_distance_matrix.txt'.format(variant_set)
  data = np.genfromtxt(f)
  np.savetxt(output_file_name, 1-data, delimiter='\t', fmt="%.09f")