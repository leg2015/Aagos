#!/bin/bash -login

# echo "in variables" TODO: make parameter list input more streamlined
# VAR1=(GENE_MOVE_PROB .001 .003 .01 .03 .1);
# VAR2=(BIT_FLIP_PROB .001 .003 .01 .03 .1);
# VAR3=(BIT_CHANGE_PROB .001 .003 .01 .03 .1); # both insert and delete
# echo "var 1 in variables script"
# echo ${VAR1[1]}

VALS_TO_TRY=( 0 .00001 .0001 .001 .003 .01 .03 .1 )