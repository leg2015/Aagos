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

# Gradient Model Runs

## N = 128

### changing bit flip rate with static environment
`bash scripts/Bash_scripts/final_launch_scripts/HpccLaunch_Mut_rates_Final.bash -VARIABLES scripts/Bash_scripts/final_launch_scripts/mut_variables_final.bash -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001  -OUTPUT_DIR Mut_Treat_Change_0_N_128 -NUM_REPLICATES 10 -GRADIENT_MODEL 1 -EMAIL [emilys email] -AAGOS_PATH [path to aagos dir]`

### f = 0.003, varying change rate

`bash scripts/Bash_scripts/final_launch_scripts/Hpcc_launch_change_env_Final.bash -VARIABLES scripts/Bash_scripts/final_launch_scripts/change_variables_final.bash -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001  -OUTPUT_DIR Change_Treat_f_.003_N_128 -NUM_REPLICATES 10 -GRADIENT_MODEL 1 -EMAIL [emilys email] -AAGOS_PATH [path to aagos dir]`

###changing bit flip rate with changing environment rate of [???]
Once an optimal environmental change rate has been discovered, run this command with said change rate
`bash scripts/Bash_scripts/final_launch_scripts/HpccLaunch_Mut_rates_Final.bash -VARIABLES scripts/Bash_scripts/final_launch_scripts/mut_variables_final.bash -GENE_MOVE_PROB .003 -BIT_INS_PROB .001 -BIT_DEL_PROB .001  -OUTPUT_DIR Mut_Treat_Change_0_N_128 -NUM_REPLICATES 10 -GRADIENT_MODEL 1 -CHANGE_RATE [???] -EMAIL [emilys email] -AAGOS_PATH [path to aagos dir]`