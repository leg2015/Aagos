
# coding: utf-8

import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
import glob
import sys
import argparse as argp
# Parses command line argument for filepath of data to clean
parser = argp.ArgumentParser(description='Clean and aggregate Aagos data.')
parser.add_argument("-f", type=str, required=True, help="filepath to where aagos data is stored. Should be the path into the dir where all run dirs are stored")
parser.add_argument("-n", type=int, required=True, help="number of replicates for this run of Aagos")
parser.add_argument("-glob", type=str, help=" what should the glob pattern for these files look like. Defaults to mutation parsing", default='/m_*/*')
args = parser.parse_args()
filepath = args.f
num_replicates = args.n
if filepath is '':
    sys.exit("No filepath was provided! Please provide filepath to Aagos data")
# script should look through fitness, representative_org and statistics file
num_files = 1 # number of files that this data script should clean out
# get path for all data files
glob_path = filepath + args.glob
files = glob.glob(glob_path)
if len(files) is 0:
    err_msg = "ERROR: files not found. Either an error with glob path: " + glob_path + " or error with data files"
    sys.exit(err_msg)
# stores all the data in this vector
dataframes_stats = []
for f in files:
        if len(files) < num_replicates:
            error_msg = "ERROR: incomplete data. A replicate from run " + f + " is missing!"
            sys.exit(error_msg)
        # get each replicate
        replicate = f.split('/')[-1]
        curr = f.split('/')[-2]
        mut_rates = curr.split('_')
        currdata = glob.glob(f + '/*.csv')
        curr_dataframes = []
        # for each file in replicate, grab the data
        for c in currdata:
            if('representative' in c): # ignore the snapshot file because not the data we're interested in right now
                curr_dataframes.append(pd.read_csv(c, index_col="update"))
            else:
                continue
        # Error check from previous issue I was having, make sure every file is present, otherwise will throw an error
        if len(curr_dataframes) < num_files:
            num_missing = (num_files - len(curr_dataframes)) 
            error_msg = "there are missing files in the data! " + str(num_missing) + " file[s] are missing from directory " + f
            sys.exit(error_msg)
            # continue
        merged = pd.concat(curr_dataframes, axis=1)
        merged["replicate"] = replicate
        for i in range(0, len(mut_rates), 2):
            merged[mut_rates[i]] = mut_rates[i+1]
        dataframes_stats.append(merged)
final = pd.concat(dataframes_stats, axis=0)
final.to_csv(filepath + '/CleanedDataRep.csv')
print("data successfully saved to ", filepath + '/CleanedDataRep.csv')
