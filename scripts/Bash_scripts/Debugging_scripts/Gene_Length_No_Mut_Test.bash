#!/bin/bash -login


INDEX=960

while [[ $INDEX -lt 1025 ]]; do
    JNDEX=0
    
    >&2 echo "---------GENE LENGTH: $INDEX ---------"
    while [[ JNDEX -lt 3 ]]; do
        time ./Aagos -GENE_MOVE_PROB 0.0 -BIT_FLIP_PROB 0.0 -BIT_INS_PROB 0.0 -BIT_DEL_PROB 0.0 -MAX_GENS 1000 -NUM_BITS $INDEX > out.txt 2>> LengthTestResults.txt
        let JNDEX=JNDEX+1
    done
    let INDEX=INDEX+1
done