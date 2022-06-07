#!/usr/bin/env python3

## Function to perform KNN on all input files and return a dictionary.
## Will return a dictionary with the 'snp_set' as key and a list of two dictionaries as values.
    ## First dictionary is the proportion of neighbours in each of the 9 populations.
    ## Second dictionary is the proportion of KNN in the correct population.
def perform_KNN(inputs, n_ind_per_pop, k):
    predictions={}
    for file in inputs:
        snp_set=re.sub(r'_[a-zA-Z_0-9.]*$','', file.split("/")[-1]) ## Infer snp_set form file name
        if snp_set == "twelve":
          snp_set = "twelve_forty"
            # ac=file.split(".")[-2][-1] ## Infer AC from file name
            # snp_set=snp_set+"_ac"+ac
        predictions[snp_set]=k_nearest_neighbour(file, n_ind_per_pop, k)
    return (predictions) 

## KNN implementation. Will return a list of 2 dictionaries. 
def k_nearest_neighbour(in_file, n_ind_per_pop, k=5):
    ## Read the input file into a numpy array
    dm = np.genfromtxt(in_file)
    
    ## Return a list of dictionaries. Shared keys. 
    ## First dictionary is the proportion of neighbours in each of the 9 populations.
    ## Second dictionary is the proportion of KNN in the correct population.
    predictions = [{},{}]
    for ind in range(len(dm)):
        ## argsort returns an array with order of indexes that would result in a sorted array. 
        ## Nearest neighbours are the indexes of the 2-(k+1)th (self-distance is always 0 and hence first). 
        population_definitions = np.repeat(["Pop0", "Pop1","Pop2","Pop3","Pop4","Pop5","Pop6","Pop7","Pop8"], n_ind_per_pop)
        nearest_neighbours = np.argsort(dm[ind])[1:k+1]
#         print(np.argsort(dm[ind])[1:k+1])
#         print(dm[ind][nearest_neighbours])
        nearest_neighbour_populations = population_definitions[nearest_neighbours]
        predictions[0]["ind"+str(ind)] = list(proportions_per_population(nearest_neighbour_populations))
        predictions[1]["ind"+str(ind)] = sum(nearest_neighbour_populations==population_definitions[ind])/k
    return (predictions)

def proportions_per_population(arr):
    proportions = []
    for i in range(9):
        proportions.append(sum(arr=="Pop"+str(i))/len(arr))
    return (proportions)

def save_results(predictions, n_ind_per_pop, m, k):
    for snp_set, (all_pops, correct_pop) in predictions.items():
        with open('KNN_K{}_{}_m{}.txt'.format(k, snp_set, m), 'w') as f:
            print("#individual", "four_mN", "variant_set", "correct_guess_proportion", *["proportion_pop"+str(i) for i in range(9)], sep="\t", file=f)
            for idx in range(n_ind_per_pop*9):
                ind="ind"+str(idx)
                print (ind, m, snp_set, correct_pop[ind], *all_pops[ind], sep="\t", file=f)

##### MAIN #####
import sys
import numpy as np
import os
import re
import matplotlib.pyplot as plt

args=sys.argv

if "-h" in args or "--help" in args or len(sys.argv) == 1:
  print("""
  A script to carry out KNN classification from a set of input distance matrices.
  All arguments are positional. Multiple distance matrices can be provided one after another.
  
  Usage:
  distance_to_KNN.py m n_ind_per_pop snp_set k input_matrix [input_matrix2 ..]
  
  Options:
    -h,--help   Print this text and exit.
  """)
  sys.exit(0)

## Read cli options.
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

snp_set=args[3]

try:
  k=int(args[4])
except ValueError:
  print("""
  ERROR: The second argument should be an integer equal to the number of sampled individuals in each of the 9 simulated populations.
  Execution halted.
  """)

input_files=args[5:]

predictions=perform_KNN(input_files, n_ind_per_pop, k)
save_results(predictions, n_ind_per_pop, m, k)