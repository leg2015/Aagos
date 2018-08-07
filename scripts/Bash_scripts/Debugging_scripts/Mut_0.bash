
#!/bin/bash -login
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -N Aagos_Mut_$START_$END
#PBS -M gillespl@southwestern.edu
#PBS -e ./$OUTPUT_DIR/error/Run_$START-$END
#PBS -o ./$OUTPUT_DIR/output/Run_$START-$END

CURR_PARAMS=""
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
-GENE_MOVE_PROB)
GENE_MOVE_PROB="$2"
# CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-SEED)
SEED="$2"
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-BIT_FLIP_PROB)
BIT_FLIP_PROB="$2"
# CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

# -BIT_INS_PROB)
# BIT_INS_PROB="$2"
# CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
# shift # shift past curr argument
# shift # shift past curr value
# ;; # indicates end of case

# -BIT_DEL_PROB)
# BIT_DEL_PROB="$2"
# CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
# shift # shift past curr argument
# shift # shift past curr value
# ;; # indicates end of case


-MAX_GENS)
MAX_GENS="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-BIT_CHANGE_PROB)
BIT_CHANGE_PROB="$2"
shift
shift
;;

-PRINT_INTERVAL)
PRINT_INTERVAL="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 " # add to list of params to include in Aagos run
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-POP_SIZE)
POP_SIZE="$2"
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

-OUTPUT_DIR) # different from DATA_FILEPATH of aagos world. that is the full path to where the data files should be stored. this val is what to name the home dir
OUTPUT_DIR="$2"
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-NUM_REPLICATES)
NUM_REPLICATES="$2"
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

    *)    # unknown option
    echo "ERROR: unknown command line argument: "$1" not a recognized command"
    exit 1
    ;;
esac
done

cd /mnt/scratch/f0004516/Aagos # make sure in Aagos directory
mkdir -vp "./$OUTPUT_DIR/"
# mkdir -vp ./$OUTPUT_DIR/m_$GENE_MOVE_PROB_f_$BIT_FLIP_PROB_c_$BIT_CHANGE_PROB/ # create dir for this replicate 
            
    let INDEX=0
        # does each replicate
        while [[ $INDEX -lt $NUM_REPLICATES ]]
        do
        # for some reason 1st time thru not being made in time or correctly
        mkdir -vp "./$OUTPUT_DIR/m_${GENE_MOVE_PROB}_f_${BIT_FLIP_PROB}_c_${BIT_CHANGE_PROB}/${INDEX}" # create dir for this replicate 
            
            FILE_PATH=./$OUTPUT_DIR/m_${GENE_MOVE_PROB}_f_${BIT_FLIP_PROB}_c_${BIT_CHANGE_PROB}/${INDEX}/ # save filepath of new dir to store dat files there
            echo $FILE_PATH
            # seed based on mutation rate and replicate. Unique seed for each run
            #  SEED=$((${GENE_MOVE_PROB}*1000+${BIT_FLIP_PROB}*100+${BIT_CHANGE_PROB}*10+${INDEX}+1))
            # SEED=1241
            # >&2 echo "\$START" # echoes mutation rate for this run to std err
            >&2 echo "Seed is: $SEED"  # echoes seed to std err
            # time run for performance purposes, pipes to std err automatically
            time ./Aagos -GENE_MOVE_PROB $GENE_MOVE_PROB -BIT_FLIP_PROB $BIT_FLIP_PROB -BIT_INS_PROB $BIT_CHANGE_PROB -BIT_DEL_PROB $BIT_CHANGE_PROB -DATA_FILEPATH $FILE_PATH -SEED $SEED $CURR_PARAMS 
            let INDEX=INDEX+1
        done
