# Smaller N runs
## N = 64
###varying bit flip mutation rates, static environment
`bash scripts/Bash_scripts/HpccLaunch_Mut_Rates_Final.bash -NUM_BITS 64 -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -NUM_REPLICATES 100 -OUTPUT_DIR Mut_Treat_Change_0_N_64 -VARIABLES scripts/Bash_scripts/mut_variables_final.bash`

### varying bit flip mutation rates, change rate of 50

`bash scripts/Bash_scripts/HpccLaunch_Mut_Rates_Final.bash -NUM_BITS 64 -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -CHANGE_RATE 50 -NUM_REPLICATES 100 -OUTPUT_DIR Mut_Treat_Change_50_N_64 -VARIABLES scripts/Bash_scripts/mut_variables_final.bash`

### f = 0.003, varying change rate

`bash scripts/Bash_scripts/Hpcc_launch_change_env.bash -NUM_BITS 64 -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -BIT_FLIP_PROB .003 -NUM_REPLICATES 100 -OUTPUT_DIR Change_Treat_f_.003_N_64 -VARIABLES scripts/Bash_scripts/change_variables.bash`


## N = 32
###varying bit flip mutation rates, static environment
`bash scripts/Bash_scripts/HpccLaunch_Mut_Rates_Final.bash -NUM_BITS 32 -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -NUM_REPLICATES 100 -OUTPUT_DIR Mut_Treat_Change_0_N_32 -VARIABLES scripts/Bash_scripts/mut_variables_final.bash`

### varying bit flip mutation rates, change rate of 50

`bash scripts/Bash_scripts/HpccLaunch_Mut_Rates_Final.bash -NUM_BITS 32 -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -CHANGE_RATE 50 -NUM_REPLICATES 100 -OUTPUT_DIR Mut_Treat_Change_50_N_32 -VARIABLES scripts/Bash_scripts/mut_variables_final.bash`

### varying bit flip mutation rates, change rate of 2

`bash scripts/Bash_scripts/HpccLaunch_Mut_Rates_Final.bash -NUM_BITS 32 -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -CHANGE_RATE 2 -NUM_REPLICATES 100 -OUTPUT_DIR Mut_Treat_Change_2_N_32 -VARIABLES scripts/Bash_scripts/mut_variables_final.bash`

### f = 0.003, varying change rate

`bash scripts/Bash_scripts/Hpcc_launch_change_env.bash -NUM_BITS 32 -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -BIT_FLIP_PROB .003 -NUM_REPLICATES 100 -OUTPUT_DIR Change_Treat_f_.003_N_32 -VARIABLES scripts/Bash_scripts/change_variables.bash`

## N = 16
###varying bit flip mutation rates, static environment
`bash scripts/Bash_scripts/HpccLaunch_Mut_Rates_Final.bash -VARIABLES scripts/Bash_scripts/mut_variables_final.bash -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -NUM_BITS 16 -OUTPUT_DIR Mut_Treat_Change_0_N_16 -NUM_REPLICATES 100`

### varying bit flip mutation rates, change rate of 50

`bash scripts/Bash_scripts/HpccLaunch_Mut_Rates_Final.bash -NUM_BITS 16 -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -CHANGE_RATE 50 -NUM_REPLICATES 100 -OUTPUT_DIR Mut_Treat_Change_50_N_16 -VARIABLES scripts/Bash_scripts/mut_variables_final.bash`

### varying bit flip mutation rates, change rate of 2

`bash scripts/Bash_scripts/HpccLaunch_Mut_Rates_Final.bash -NUM_BITS 16 -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -CHANGE_RATE 2 -NUM_REPLICATES 100 -OUTPUT_DIR Mut_Treat_Change_2_N_16 -VARIABLES scripts/Bash_scripts/mut_variables_final.bash`

### f = 0.003, varying change rate

`bash scripts/Bash_scripts/Hpcc_launch_change_env.bash -VARIABLES scripts/Bash_scripts/change_variables.bash -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001 -NUM_BITS 16 -OUTPUT_DIR Change_Treat_f_.003_N_16 -NUM_REPLICATES 10`