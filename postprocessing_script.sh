#!/bin/bash

# Makes postprocessing figures

runname=$1
echo "Runname is "$runname

if test -f "/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/execute_model.py"; 
then 
echo "WE ARE WORKING IN MAC MODE"
echo "Making deformation data comparisons"
python /Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/compare_with_deformation_data.py $runname
python /Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/compare_with_deformation_data.py $runname deephead=True
echo "Plotting aquifer partitioning"
python /Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/plot_partitioning.py $runname
echo "Plotting interbed partitioning"
python /Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/plot_partitioning_interbeds.py $runname
fi

if test -f "/home/mlees/Land_Subsidence/Local_Scale/compaction_model/execute_model.py";
then 
echo "WE ARE WORKING IN LINIUX MODE"
echo "Making deformation data comparisons"
python /home/mlees/Land_Subsidence/Local_Scale/compaction_model/compare_with_deformation_data.py $runname
python /home/mlees/Land_Subsidence/Local_Scale/compaction_model/compare_with_deformation_data.py $runname deephead=True
echo "Plotting aquifer partitioning"
python /home/mlees/Land_Subsidence/Local_Scale/compaction_model/plot_partitioning.py $runname
echo "Plotting interbed partitioning"
python /home/mlees/Land_Subsidence/Local_Scale/compaction_model/plot_partitioning_interbeds.py $runname
fi

if test -f "/home/mlees/compaction_model/execute_model.py";
then
echo "WE ARE WORKING IN Knightblade MODE"
echo "Making deformation data comparisons"
python /home/mlees/compaction_model/compare_with_deformation_data.py $runname
python /home/mlees/compaction_model/compare_with_deformation_data.py $runname deephead=True
echo "Plotting aquifer partitioning"
python /home/mlees/compaction_model/plot_partitioning.py $runname
echo "Plotting interbed partitioning"
python /home/mlees/compaction_model/plot_partitioning_interbeds.py $runname
fi


