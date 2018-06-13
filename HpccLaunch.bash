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
CURR_PARAMS=" $CURR_PARAMS $key $2 "
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
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-MAX_GENS)
MAX_GENS="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
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
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-TOURNAMENT_SIZE)
TOURNAMENT_SIZE="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-NUM_GENES)
NUM_GENES="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-NUM_BITS)
NUM_BITS="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-GENE_SIZE)
GENE_SIZE="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-MAX_SIZE)
MAX_SIZE="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-MIN_SIZE)
MIN_SIZE="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-GENE_MOVE_PROB)
GENE_MOVE_PROB="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-BIT_FLIP_PROB)
BIT_FLIP_PROB="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-BIT_INS_PROB)
BIT_INS_PROB="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-BIT_DEL_PROB)
BIT_DEL_PROB="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-PRINT_INTERVAL)
PRINT_INTERVAL="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-STATISTICS_INTERVAL)
STATISTICS_INTERVAL="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-SNAPSHOT_INTERVAL)
SNAPSHOT_INTERVAL="$2"
CURR_PARAMS=" $CURR_PARAMS $key $2 "
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-NUM_REPLICATES)
NUM_REPLICATES="$2"
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-NUM_JOBS)
NUM_JOBS="$2"
shift # shift past curr argument
shift # shift past curr value
;; # indicates end of case

-NUM_CONDITIONS)
NUM_CONDITIONS="$2"
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

# if [[ $VARIABLES -eq "NULL" ]]; then # TODO: this is complaining
#    echo "ERROR: must have a file with changing variables as a command line argument"
#    exit 1
# fi

if [[ $NUM_REPLICATES -eq "NULL" ]]; then # TODO: this is complaining
   echo "no replicate given, defaulting to 10"
   let NUM_REPLICATES=10
fi

if [[ $SEED -eq "NULL" ]]; then # TODO: this is complaining
   echo "no Seed given, defaulting to 1"
   let SEED=1
fi


# now have access to the changing parameters!!!
source ${VARIABLES}

NUM_VALS=${#VALS_TO_TRY[@]}
NUM_VARS=${#ALL_VARS[@]}
TOT_RUNS=$(($NUM_VALS ** $NUM_VARS))



#TODO: assuming bash script and aagos excecutable are in same directory, maybe add a check to confirm?
# cd ..


mkdir -vp "./$OUTPUT_DIR/scripts" 

# cd  "./$OUTPUT_DIR" # TODO: get rid of this move too!!


# mkdir -v "./scripts"

# cd "./scripts"



CURR=0

# TODO: right now loop only works for the parameters we're checking right now. Not sure how to generalize, need to figure out
# CURR_IND=()
# while [[ CURR -lt $NUM_VARS ]]; do
# CURR_IND+=(0)
#  let CURR=CURR+1
# done

# CURR_IND=(0 0 0) #TODO: this is what's inflexible currently
# echo ${CURR_IND[0]}
# echo ${CURR_IND[1]}
# echo ${CURR_IND[2]}
# echo "number of jobs "
# echo $NUM_JOBS

# -------------------- combinatorics that works ----------------------------- #
IND_0=0
IND_1=0
IND_2=0
SIZE=${#VALS_TO_TRY[@]}
COUNT=0
# echo "current directory before loop: "$PWD #aagos
while [[ $IND_0 -lt  SIZE ]]; do
    let CURR_IND[0]=$IND_0
    let CURR_IND[1]=$IND_1
    let CURR_IND[2]=$IND_2
    IND_TEMP=0
    let CURR_IND[1]=IND_TEMP
    let END_IND=$((SIZE-1))
    START=m_${VALS_TO_TRY[CURR_IND[0]]}_f_${VALS_TO_TRY[CURR_IND[1]]}_c_${VALS_TO_TRY[CURR_IND[2]]}
    END=m_${VALS_TO_TRY[CURR_IND[0]]}_f_${VALS_TO_TRY[END_IND]}_c_${VALS_TO_TRY[CURR_IND[2]]} # TODO: not sure if arithmetic is inlinable
# echo "current directory: "$PWD" count:"$COUNT
let COUNT=COUNT+1
INDEX=0

# echo "befire"$PWD

# TODO: move bash calls instead, don't have to be in the directory to run as long as filepath call is correct!!
# if we run the sub files from output dir then don't have to do weird cd s

# echo "num rep "$NUM_REPLICATES
# TODO: assuming in output dir
# cat << EOF > "./$OUTPUT_DIR/scripts/Run_$START-$END.sub" #TODO: fix other cat file!!
cat << EOF > "./$OUTPUT_DIR/scripts/Run_$START-$END.qsub" #TODO: fix other cat file!!
#!/bin/bash -login
#PBS -l walltime=00:04:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -N Aagos_Mut_$START_$END
#PBS -M gillespl@southwestern.edu
#PBS -j oe
#PBS -o ./OE_Run_$START-$END

cd /mnt/scratch/f0004516/Aagos

source  variables.bash


let IND_TEMP=$IND_TEMP
    # echo "first cat "\${CURR_IND[1]}
    # echo "first val"\${VALS_TO_TRY[CURR_IND[1]]}
    while [[ \$IND_TEMP -lt $SIZE ]]; do   # TODO: will the 0 thats ind temp break everything?
# echo \$(((${CURR_IND[0]}+1)*100+(\${CURR_IND[1]}+1)*10+(\${CURR_IND[2]}+1))) 
       
        START=m_${VALS_TO_TRY[CURR_IND[0]]}_f_\${VALS_TO_TRY[CURR_IND[1]]}_c_\${VALS_TO_TRY[CURR_IND[2]]}   
        let INDEX=0
        while [[ \$INDEX -lt $NUM_REPLICATES ]]; do
        mkdir -vp "./$OUTPUT_DIR/\$START/\$INDEX" #TODO: change
        FILE_PATH=./$OUTPUT_DIR/\$START/\$INDEX/ # TODO: change

# TODO: make sure file path given to aagos works correctly

let SEED=\$(((${CURR_IND[0]}+1)*1000+(\${CURR_IND[1]}+1)*100+(\${CURR_IND[2]}+1)*10+\$INDEX+1)) 
# echo \$SEED


# TODO: im not sure the params below are correct, may need to escape them....
 ./Aagos -GENE_MOVE_PROB ${VALS_TO_TRY[CURR_IND[0]]} -BIT_FLIP_PROB \${VALS_TO_TRY[CURR_IND[1]]} -BIT_INS_PROB \${VALS_TO_TRY[CURR_IND[2]]} -BIT_DEL_PROB \${VALS_TO_TRY[CURR_IND[2]]} -DATA_FILEPATH \$FILE_PATH -SEED \$SEED $CURR_PARAMS # TODO: assumes using the rest of vals as default
        
        let INDEX=INDEX+1
        done
        let IND_TEMP=IND_TEMP+1
        let CURR_IND[1]=\$IND_TEMP
        
    done


EOF
# let SEED=$(((\${CURR_IND[0]}+1)*1000+(\${CURR_IND[1]}+1)*100+(\${CURR_IND[2]}+1)*10+\$INDEX+1)) 

# echo "after first cat "$PWD
#in AAGOS

#  echo "sub "./Run_$START-$END.sub"" # TODO: change echo
qsub ./$OUTPUT_DIR/scripts/Run_$START-$END.qsub
    IND_2=1
        let CURR_IND[0]=$IND_0
    let CURR_IND[1]=$IND_1
    let CURR_IND[2]=$IND_2
    let INDEX=0
    while [[ $IND_1 -lt SIZE ]]; do
                START=m_${VALS_TO_TRY[CURR_IND[0]]}_f_${VALS_TO_TRY[CURR_IND[1]]}_c_${VALS_TO_TRY[CURR_IND[2]]}
                END=m_${VALS_TO_TRY[CURR_IND[0]]}_f_${VALS_TO_TRY[CURR_IND[1]]}_c_${VALS_TO_TRY[$END_IND]]} 
            # echo "current dir in 2nd while loop: "$PWD" count:"$COUNT
            let COUNT=COUNT+1

            memes="memes"
# echo "cur vals: $CURR_PARAMS"
# cat << EOF > "./$OUTPUT_DIR/scripts/Run_$START-$END.sub"
cat << EOF > "./$OUTPUT_DIR/scripts/Run_$START-$END.qsub" #TODO: check that script can run correctly
#!/bin/bash -login
#PBS -l walltime=00:04:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -N Aagos_Mut_$START_$END
#PBS -M gillespl@southwestern.edu
#PBS -j oe
#PBS -o ./OE_Run_$START-$END

cd /mnt/scratch/f0004516/Aagos
source  variables.bash

        let IND_2=$IND_2
        while [[ \$IND_2 -lt $SIZE ]]; do #TODO: will the ind 2 break?
            let CURR_IND[0]=$IND_0
            let CURR_IND[1]=$IND_1
            let CURR_IND[2]=\$IND_2
            START=m_\${VALS_TO_TRY[CURR_IND[0]]}_f_\${VALS_TO_TRY[CURR_IND[1]]}_c_\${VALS_TO_TRY[CURR_IND[2]]}   

        #  echo \$(((\${CURR_IND[0]}+1)*100+(\${CURR_IND[1]}+1)*10+(\${CURR_IND[2]}+1))) 
       
        let INDEX=0

        while [[ \$INDEX -lt $NUM_REPLICATES ]]; do

 
        mkdir -vp "./$OUTPUT_DIR/\$START/\$INDEX"
        FILE_PATH=./$OUTPUT_DIR/\$START/\$INDEX/ # TODO: see if this works
        let SEED=\$(((\${CURR_IND[0]}+1)*1000+(\${CURR_IND[1]}+1)*100+(\${CURR_IND[2]}+1)*10+\$INDEX+1)) 
        # echo \$SEED
         ./Aagos -GENE_MOVE_PROB \${VALS_TO_TRY[CURR_IND[0]]} -BIT_FLIP_PROB \${VALS_TO_TRY[CURR_IND[1]]} -BIT_INS_PROB \${VALS_TO_TRY[CURR_IND[2]]} -BIT_DEL_PROB \${VALS_TO_TRY[CURR_IND[2]]} -DATA_FILEPATH \$FILE_PATH -SEED \$SEED $CURR_PARAMS # TODO: assumes using the rest of vals as default
 
               let INDEX=INDEX+1
        done


    let IND_2=IND_2+1  
    done

EOF
    #    echo "sub "./Run_$START-$END.sub""# TODO: change echo   

    qsub ./$OUTPUT_DIR/scripts/Run_$START-$END.qsub

    

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
let COUNT=COUNT+1
# -------------------- combinatorics that works ----------------------------- #

# cat << EOF > "Run_$START-$END.sub"

# #!/bin/bash -login

# #PBS -l walltime=00:04:00
# #PBS -l nodes=1:ppn=1
# #PBS -l mem=2gb
# #PBS -N Aagos_Mut_$START_$END
# #PBS -M gillespl@southwestern.edu
# #PBS -j oe
# #PBS -o ./OE_Run_$START-$END

# SEED=$SEED

# cd .. # move to $OUTPUT_DIR/

# # cd .. # move to Aagos/
#  # now just need to actually run this version of Aagos 4-5 times

# START=m_${VALS_TO_TRY[CURR_IND[0]]}_f_${VALS_TO_TRY[CURR_IND[1]]}_c_${VALS_TO_TRY[CURR_IND[2]]}

# CURR=0
#  while [[]]; do # TODO: figure out while loop, updating values correctly
# mkdir -v "./$START"
# IN_CURR=0

# while [[ $IN_CURR -lt $NUM_REPLICATES ]]; do
# mkdir -v "./$START_REP_$IN_CURR"
# cd "./$START_REP_$IN_CURR"
# FILE_PATH=pwd # TODO: see if this works
# echo $FILE_PATH
# cd .. 
# cd ..
# cd ..
# pwd
#  # ./Aagos -GENE_MOVE_PROB ${VALS_TO_TRY[CURR_IND[0]]} -BIT_FLIP_PROB ${VALS_TO_TRY[CURR_IND[1]]} -BIT_INS_PROB ${VALS_TO_TRY[CURR_IND[2]]} -BIT_DEL_PROB ${VALS_TO_TRY[CURR_IND[2]]} -OUTPUT_DIR $FILE_PATH - $START -SEED $SEED
# done


# let SEED=SEED+1 # TODO: see if seed increments correctly
# done # TODO: include 10 replicates for each run

# EOF
# # make sure there is NO OTHER character with the EOF above!
# let CURR=CURR+1
# # echo "seed after"
# # echo $SEED
# let SEED=SEED+1 # TODO: this seed update needs to change
# done

# echo "final seed is "
# echo $SEED


# steps: 
# 1. create bash file that in turn
#     2. creates a bunch of *.qsub files
#     3. each .qsub file looks like the example file given online
#     4. that's where the make call and everything will be made


