#!/bin/bash -login

# parameter parsing code derived from this example:
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash 
CURR_PARAMS=""
while [[ $# -gt 0 ]]
do
key="$1"

# syntax to access vars is ${var_name}
case $key in

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

-)BIT_INS_PROB
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
if [[ $VARIABLES -eq "NULL" ]]; then 
   echo "ERROR: must have a file with changing variables as a command line argument"
   exit 1
fi
# must have a val for num replicates
if [[ $NUM_REPLICATES -eq "NULL" ]]; then 
   echo "no replicate given, defaulting to 10"
   let NUM_REPLICATES=10
fi

# Gives access to changing params
source ${VARIABLES}
NUM_ENV=${#CHANGE_TO_TRY[@]}


#TODO: assuming bash script and aagos excecutable are in same directory, maybe add a check to confirm?
mkdir -vp "./$OUTPUT_DIR/scripts" 
mkdir -vp "./$OUTPUT_DIR/output"
mkdir -vp "./$OUTPUT_DIR/error"
# Sets up vars that keep track of what mutation to look at currently
IND=0
# This while loop actually creates the scripts that run on the hpcc.
# Creates 30 scripts which each launch 10 replicates of either 4 or 5 mutation rate combos
while [[ $IND -lt  NUM_ENV ]]; do
    START=change_${CHANGE_TO_TRY[$IND]}
    INDEX=0
    # This cat creates new script in script dir. Assumes in Aagos dir.
cat << EOF > "./$OUTPUT_DIR/scripts/Run_$START.qsub"
#!/bin/bash -login
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -N Aagos_Change_$START
#PBS -M gillespl@southwestern.edu
#PBS -e ./$OUTPUT_DIR/error/Run_$START
#PBS -o ./$OUTPUT_DIR/output/Run_$START

cd /mnt/scratch/f0004516/Aagos # make sure in Aagos directory
source  $VARIABLES # variables must be in Aagos dir to work
let INDEX=0
while [[ \$INDEX -lt $NUM_REPLICATES ]]
    do
        START=change_${CHANGE_TO_TRY[$IND]}
        # does each replicate
        mkdir -vp "./$OUTPUT_DIR/$START/\$INDEX" # create dir for this replicate 
        FILE_PATH=./$OUTPUT_DIR/$START/\$INDEX/ # save filepath of new dir to store dat files there
        # seed based on mutation rate and replicate. Unique seed for each run
        let SEED=\$((${CHANGE_TO_TRY[$IND]}*100+\$INDEX+1))
        >&2 echo "\$START" # echoes mutation rate for this run to std err
        >&2 echo "\$SEED"  # echoes seed to std err
        # time run for performance purposes, pipes to std err automatically
        time ./Aagos -DATA_FILEPATH \$FILE_PATH -SEED \$SEED $CURR_PARAMS -CHANGE_RATE ${CHANGE_TO_TRY[$IND]}
        let INDEX=INDEX+1
    done
EOF

    # should still be in Aagos here
    # launch scripts to hpcc scheduler to be run
    qsub "./$OUTPUT_DIR/scripts/Run_$START.qsub"
    let IND=IND+1
done