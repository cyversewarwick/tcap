#!/bin/bash
set -e

#average reps
python3 /scripts/expression_filter.py --Input $1 --AverageReps

#run TCAP, passing all the other arguments
python3 /scripts/tcap.py --Input FilteredOutput.csv "${@:2}"

#make folder for bingo/meme output
mkdir functional_analysis_inputs

#make bingo/meme output
python3 /scripts/bingomeme.py

#kick out the averaged reps expression file
rm FilteredOutput.csv