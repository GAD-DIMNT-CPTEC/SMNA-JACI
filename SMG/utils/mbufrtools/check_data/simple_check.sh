#!/bin/bash
source ./dirconf.sh
clear
echo "------------------------------------------------"
echo " simple_check                                   " 
echo "------------------------------------------------"
echo "INPUT_DIR = "$INPUT_DIR
echo "OUTPUT_DIR = "$OUTPUT_DIR
source ./get_parameters.sh 

#INPUT=${INPUT_DIR}/$yy0/$mm0/$dd0
#echo $INPUT
#ls -ltr $INPUT

LIST=${INPUT_DIR}/$yy0/$mm0/$dd0/*${zz0}z*
../bin/bufrcontent $LIST

