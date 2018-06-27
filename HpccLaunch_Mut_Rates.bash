#!/bin/bash -login

# parameter parsing code derived from this example:
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash 
CURR_PARAMS=""
while [[ $# -gt 0 ]]
do
key="$1"

# syntax to access vars is ${var_name}
case $key in
-CHANGE_RATE)
CHANGE_RATE="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case


-OUTPUT_DIR) # different from DATA_FILEPATH of aagos world. that is the full path to where the data files should be stored. this val is what to name the home dir
OUTPUT_DIR="$2"
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-POP_SIZE)
POP_SIZE="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-MAX_GENS)
MAX_GENS="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-SEED)
SEED="$2"
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case


-ELITE_COUNT)
ELITE_COUNT="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-TOURNAMENT_SIZE)
TOURNAMENT_SIZE="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-NUM_GENES)
NUM_GENES="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-NUM_BITS)
NUM_BITS="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-GENE_SIZE)
GENE_SIZE="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-MAX_SIZE)
MAX_SIZE="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-MIN_SIZE)
MIN_SIZE="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-GENE_MOVE_PROB)
GENE_MOVE_PROB="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-BIT_FLIP_PROB)
BIT_FLIP_PROB="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-BIT_INS_PROB)
BIT_INS_PROB="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-BIT_DEL_PROB)
BIT_DEL_PROB="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-PRINT_INTERVAL)
PRINT_INTERVAL="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-STATISTICS_INTERVAL)
STATISTICS_INTERVAL="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-SNAPSHOT_INTERVAL)
SNAPSHOT_INTERVAL="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-NUM_REPLICATES)
NUM_REPLICATES="$2"
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-VARIABLES)
VARIABLES="$2"
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case



    *)    # unknown option
    echo "ERROR: unknown command line argument: "$1" not a recognized command"
    exit 1
    ;;
esac
done
# must have a var file for script to run properly
# if [[ $VARIABLES -eq "NULL" ]]; then 
#    echo "ERROR: must have a file with changing variables as a command line argument"
#    exit 1
# fi
# must have a val for num replicates
if [[ $NUM_REPLICATES -eq "NULL" ]]; then 
   echo "no replicate given, defaulting to 10"
   let NUM_REPLICATES=10
fi

# Gives access to changing params
source ${VARIABLES}
# stores vals from var file
NUM_VALS=${#VALS_TO_TRY[@]}
NUM_VARS=${#ALL_VARS[@]}
TOT_RUNS=$(($NUM_VALS ** $NUM_VARS))



#TODO: assuming bash script and aagos excecutable are in same directory, maybe add a check to confirm?
mkdir -vp "./$OUTPUT_DIR/scripts" 
mkdir -vp "./$OUTPUT_DIR/output"
mkdir -vp "./$OUTPUT_DIR/error"
# Sets up vars that keep track of what mutation to look at currently
IND_0=0
IND_1=0
IND_2=0
SIZE=${#VALS_TO_TRY[@]}
COUNT=0
# This while loop actually creates the scripts that run on the hpcc.
# Creates 30 scripts which each launch 10 replicates of either 4 or 5 mutation rate combos
while [[ $IND_0 -lt  SIZE ]]; do
    # reset what mutation combo we're looking at right now
    let CURR_IND[0]=$IND_0
    let CURR_IND[1]=$IND_1
    let CURR_IND[2]=$IND_2
    IND_TEMP=0
    let CURR_IND[1]=IND_TEMP # first loop changes just second mutation rate
    let END_IND=$((SIZE-1))
    # gets what mutation rates are being explored in this script
    START=m_${VALS_TO_TRY[CURR_IND[0]]}_f_${VALS_TO_TRY[CURR_IND[1]]}_c_${VALS_TO_TRY[CURR_IND[2]]}
    END=m_${VALS_TO_TRY[CURR_IND[0]]}_f_${VALS_TO_TRY[END_IND]}_c_${VALS_TO_TRY[CURR_IND[2]]} 
    let COUNT=COUNT+1
    INDEX=0
    # This cat creates new script in script dir. Assumes in Aagos dir.
cat << EOF > "./$OUTPUT_DIR/scripts/Run_$START-$END.qsub"
#!/bin/bash -login
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -N Aagos_Mut_$START_$END
#PBS -M gillespl@southwestern.edu
#PBS -e ./$OUTPUT_DIR/error/Run_$START-$END
#PBS -o ./$OUTPUT_DIR/output/Run_$START-$END

cd /mnt/scratch/f0004516/Aagos # make sure in Aagos directory
source  $VARIABLES # variables must be in Aagos dir to work
let IND_TEMP=$IND_TEMP # save var from variables file

    while [[ \$IND_TEMP -lt $SIZE ]]
    do
        START=m_${VALS_TO_TRY[CURR_IND[0]]}_f_\${VALS_TO_TRY[CURR_IND[1]]}_c_\${VALS_TO_TRY[CURR_IND[2]]}   
        let INDEX=0
        # does each replicate
        while [[ \$INDEX -lt $NUM_REPLICATES ]]
        do
            mkdir -vp "./$OUTPUT_DIR/\$START/\$INDEX" # create dir for this replicate 
            FILE_PATH=./$OUTPUT_DIR/\$START/\$INDEX/ # save filepath of new dir to store dat files there
            # seed based on mutation rate and replicate. Unique seed for each run
            let SEED=\$(((${CURR_IND[0]}+1)*1000+(\${CURR_IND[1]}+1)*100+(\${CURR_IND[2]}+1)*10+\$INDEX+1))
            >&2 echo "\$START" # echoes mutation rate for this run to std err
            >&2 echo "\$SEED"  # echoes seed to std err
            # time run for performance purposes, pipes to std err automatically
            time ./Aagos -GENE_MOVE_PROB ${VALS_TO_TRY[CURR_IND[0]]} -BIT_FLIP_PROB \${VALS_TO_TRY[CURR_IND[1]]} -BIT_INS_PROB \${VALS_TO_TRY[CURR_IND[2]]} -BIT_DEL_PROB \${VALS_TO_TRY[CURR_IND[2]]} -DATA_FILEPATH \$FILE_PATH -SEED \$SEED $CURR_PARAMS 
            let INDEX=INDEX+1
        done
        #updates loop variables
        let IND_TEMP=IND_TEMP+1
        let CURR_IND[1]=\$IND_TEMP
    done
EOF

    # should still be in Aagos here
    # launch scripts to hpcc scheduler to be run
    qsub "./$OUTPUT_DIR/scripts/Run_$START-$END.qsub"
    # reset second mut_rate, will always be starting at this val
    IND_2=1
    let CURR_IND[0]=$IND_0
    let CURR_IND[1]=$IND_1
    let CURR_IND[2]=$IND_2
    # reset index var
    let INDEX=0
    while [[ $IND_1 -lt SIZE ]]
    do
        START=m_${VALS_TO_TRY[CURR_IND[0]]}_f_${VALS_TO_TRY[CURR_IND[1]]}_c_${VALS_TO_TRY[CURR_IND[2]]}
        END=m_${VALS_TO_TRY[CURR_IND[0]]}_f_${VALS_TO_TRY[CURR_IND[1]]}_c_${VALS_TO_TRY[$END_IND]]} 
        let COUNT=COUNT+1
cat << EOF > "./$OUTPUT_DIR/scripts/Run_$START-$END.qsub"
#!/bin/bash -login
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -N Aagos_Mut_$START_$END
#PBS -M gillespl@southwestern.edu
#PBS -e ./$OUTPUT_DIR/error/Run_$START-$END
#PBS -o ./$OUTPUT_DIR/output/Run_$START-$END

cd /mnt/scratch/f0004516/Aagos # make sure in Aagos directory
source  $VARIABLES # variables must be in Aagos dir to work
let IND_2=$IND_2 # save var from variables file
        
while [[ \$IND_2 -lt $SIZE ]]
do 
    let CURR_IND[0]=$IND_0
    let CURR_IND[1]=$IND_1
    let CURR_IND[2]=\$IND_2
    START=m_\${VALS_TO_TRY[CURR_IND[0]]}_f_\${VALS_TO_TRY[CURR_IND[1]]}_c_\${VALS_TO_TRY[CURR_IND[2]]}   
    let INDEX=0
    # runs each replicate
    while [[ \$INDEX -lt $NUM_REPLICATES ]]
    do
        mkdir -vp "./$OUTPUT_DIR/\$START/\$INDEX" # creates dir for current replicate
        FILE_PATH=./$OUTPUT_DIR/\$START/\$INDEX/ 
        let SEED=\$(((\${CURR_IND[0]}+1)*1000+(\${CURR_IND[1]}+1)*100+(\${CURR_IND[2]}+1)*10+\$INDEX+1)) 
        >&2 echo "\$START" # pipes mut_rate to std err
        >&2 echo "\$SEED"  # pipes seed to std err
        # time aagos run for performance 
         time ./Aagos -GENE_MOVE_PROB \${VALS_TO_TRY[CURR_IND[0]]} -BIT_FLIP_PROB \${VALS_TO_TRY[CURR_IND[1]]} -BIT_INS_PROB \${VALS_TO_TRY[CURR_IND[2]]} -BIT_DEL_PROB \${VALS_TO_TRY[CURR_IND[2]]} -DATA_FILEPATH \$FILE_PATH -SEED \$SEED $CURR_PARAMS 
        let INDEX=INDEX+1
    done
    let IND_2=IND_2+1  
done
EOF
        qsub "./$OUTPUT_DIR/scripts/Run_$START-$END.qsub"   
        # update all variables
        let IND_2=1
        let IND_1=IND_1+1
        let CURR_IND[0]=$IND_0
        let CURR_IND[1]=$IND_1
        let CURR_IND[2]=$IND_2
    done
    let IND_1=0 
    let IND_2=0
    let IND_0=IND_0+1
done