#!/bin/bash

# Makes postprocessing figures

runname=$1
echo "Runname is "$runname

echo "Making deformation data comparisons"
python /Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/compare_with_deformation_data.py $runname
python /Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/compare_with_deformation_data.py $runname deephead=True

echo "Plotting aquifer partitioning"
python /Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/plot_partitioning.py $runname


echo "Plotting interbed partitioning"
python /Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/plot_partitioning_interbeds.py $runname