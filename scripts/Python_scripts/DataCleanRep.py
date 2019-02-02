
# coding: utf-8
# how to run this: run from aagos direcotyr: python [path to file]/DataCleanRep.py -f [path to group of files to clean] -n [number of replicates for that trial] -glob [the path to the data tracking files ie representative.csv etc.]
# glob path for change files is [ /c*/* ]
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
parser.add_argument("-glob", type=str, help=" what should the glob pattern for these files look like. Defaults to mutation parsing", default='/m_*/*') # default is for using with mutation tests bash scripts
parser.add_argument("-change", action="store_true", help=" if the file that's being cleaned is for changing environments") # if using changing environments bash script, then use this flag and will automatically pull changing environemnts for you

args = parser.parse_args()
filepath = args.f
filename = ''
num_replicates = args.n
if args.change is True:
    glob_path = filepath + "/c*/*"
else :
    glob_path = filepath + args.glob
if filepath is '':
    sys.exit("No filepath was provided! Please provide filepath to Aagos data")
# script should look through fitness, representative_org and statistics file
num_files = 1 # number of files that this data script should clean out
# get path for all data files

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
        # print("f is: " + str(f))
        # print("f at 0 is " + str(f.split('/')[1])) # works as long as run from Aagos directory
        filename = f.split('/')[1] # only works if python script run from Aagos directory!!
        replicate = f.split('/')[-1]
        curr = f.split('/')[-2]
        mut_rates = curr.split('_')
        currdata = glob.glob(f + '/*.csv')
        curr_dataframes = []
        # print("replicate is " + str(replicate))
        # print("curr is: " + str(curr))
        # print("mut rates is " + str(mut_rates))
        # print("curr data is: " + str(currdata))
        # for each file in replicate, grab the data
        for c in currdata:
            if('representative' in c): # ignore the snapshot file because not the data we're interested in right now
                curr_dataframes.append(pd.read_csv(c, index_col="update"))
                print("Found the representative file")
            else:
                continue
        # Error check from previous issue I was having, make sure every file is present, otherwise will throw an error
        if len(curr_dataframes) < num_files:
            num_missing = (num_files - len(curr_dataframes)) 
            print("The number of expected files is: " + str(num_files) + " and the number of files found is: " + str(len(curr_dataframes)))
            error_msg = "there are missing files in the data! " + str(num_missing) + " file[s] are missing from directory " + f
            sys.exit(error_msg)
            # continue
        merged = pd.concat(curr_dataframes, axis=1)
        merged["replicate"] = replicate
        for i in range(0, len(mut_rates), 2):
            merged[mut_rates[i]] = mut_rates[i+1]
        dataframes_stats.append(merged)
final = pd.concat(dataframes_stats, axis=0)
final.to_csv(filepath +  '/' + filename + '_CleanedDataRep.csv')
print("data successfully saved to ", filepath +  '/' + filename + '_CleanedDataRep.csv')
