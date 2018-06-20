
# coding: utf-8

# In[85]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
import argparse as argp

parser = argp.ArgumentParser(description='Clean and aggregate Aagos data.')
parser.add_argument("-f", type=str, help="filepath to where aagos data is stored. Should be the path into the dir where all run dirs are stored")
# parser.add_argument('f', metavar='filepath', type='string',  default=argparse.SUPPRESS, help='filepath to where aagos data is stored. Should be the path into the dir where all run dirs are stored')
args = parser.parse_args()
filepath = args.f
if filepath is None:
    sys.exit("No filepath was provided! Please provide filepath to Aagos data")
print("filepath is ", filepath)
print(args)


# In[71]:


num_files = 3 # number of files that this data script should clean out
files = glob.glob(filepath + '/m_*/*')
dataframes_stats = []
for f in files:
        replicate = f.split('/')[-1]
        curr = f.split('/')[-2]
        mut_rates = curr.split('_')
        currdata = glob.glob(f + '/*.csv')
        curr_dataframes = []
        for c in currdata:
            if('snapshot' in c):
                continue
            curr_dataframes.append(pd.read_csv(c, index_col="update"))
    
        if len(curr_dataframes) < num_files: # so c is left over from the previous file set! thats why its there
            num_missing = (num_files - len(curr_dataframes)) 
            error_msg = "there are missing files in the data! " + str(num_missing) + " file[s] are missing from directory " + f
            sys.exit(error_msg)
            continue
        merged = pd.concat(curr_dataframes, axis=1)
        merged["replicate"] = replicate
        for i in range(0, len(mut_rates), 2):
            merged[mut_rates[i]] = mut_rates[i+1]
        dataframes_stats.append(merged)
final = pd.concat(dataframes_stats, axis=0)
final.to_csv(filepath + '/CleanedData.csv')

