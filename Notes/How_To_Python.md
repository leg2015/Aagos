# How to run the python script

## DataCleanRep.py, DataCleanStatFit.py, and DataCleanStatFitRep.py:s
get in Aagos directory and run:
### For changing mutation rates runs: 
'python scripts/Python_scripts/[script name].py -f [filepath to directory holding all the runs youre interested in] -n 10[could change but usually this] -glob m_.003_f_.*_c_.001/*'

### For changing environment rate runs:
'python scripts/Python_scripts/[Script name].py -f [[filepath to directory holding all the runs youre interested in]] -n 10 -glob change_*/*'

## DataCleanStat.py
this one is broken currently, hold off on for now...