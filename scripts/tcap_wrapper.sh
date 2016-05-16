#!/bin/bash

#run TCAP, passing all the arguments
python3 /scripts/tcap.py "${@:1}"

#make folder for bingo/meme output
mkdir functional_analysis_inputs

#make bingo/meme output
python3 bingomeme.py