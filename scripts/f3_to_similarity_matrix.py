#!/usr/bin/env python3 

### Function to read in a log file and update the f3_matrix
def make_f3_matrix(f3_matrix, in_file):
    data = np.genfromtxt(in_file, autostrip=True, dtype=str, skip_header=13, usecols=[1,2,4])
    data = np.char.replace(data, 'ind', '').astype(np.float32)
    for row in data:
        f3_matrix[int(row[0])][int(row[1])]+=row[2]/1000
        f3_matrix[int(row[1])][int(row[0])]+=row[2]/1000

    return(f3_matrix)

import sys
import numpy as np

args=sys.argv[1:]

if "-h" in args or "--help" in args or len(sys.argv) == 1:
  print("""
  A small script to compile a similarity matrix out of a provided set of qp3Pop output files, plus the number of individuals in each population. Assumes there are 9 populations.
  All arguments are positional. Multiple log files can be provided one after another.
  
  Usage:
  f3_to_distance_matrix.py <N_per_pop> <output_file_name> <qp3Pop Logfiles>
  
  Options:
    -h,--help   Print this text and exit.
  """)
  sys.exit(0)


try:
  n_per_pop=int(args[0])
except ValueError:
  print("""
  ERROR: The first argument should be an integer equal to the number of sampled individuals in each of the 9 simulated populations.
  Execution halted.
  """)

output_file_name=args[1]

f3_logs=args[2:]

num_inds = n_per_pop*9
f3_matrix = np.zeros((num_inds,num_inds),dtype=np.float64)

for log in f3_logs:
  f3_matrix=make_f3_matrix(f3_matrix, log)


## Since the values within the matrix are the sum of f3 across all chromosomes, the matrix then needs to be divided by the number of chromosomes to get the average f3 across all chromosomes.
## After division, the diagonal needs to be reset to 1 to reflect 0 distance within individual.

f3_matrix=f3_matrix/20
np.fill_diagonal(f3_matrix, 1)

np.savetxt(output_file_name, f3_matrix, delimiter='\t', fmt="%.07f")

# for log in f3_logs:
#   for line in open(log, "r"):
#     fields=line.strip().split()
#     if fields[0] != "result:":
#       continue
#     indA=int(fields[1].lstrip("ind"))
#     indB=int(fields[2].lstrip("ind"))
#     ## 'outgroupmode: YES' sets the f3 denominator to 0.001, which we need to reverse by dividing the value given by 1000
#     f3_value=float(fields[4]) / 1000
#     print(indA, indB, '{:.9f}'.format(f3_value))