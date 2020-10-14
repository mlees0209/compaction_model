#!/bin/bash

# Input arguments are: 1) a master file of format::
#
# runfamilyname
# dummy_paramfilename
# variable to change
# run names
# variable values (separated by |)

masterfile=$1
echo "Masterfile found to be "$masterfile

# Read in the inputs
dummyfilename=$(head -n 2 $masterfile | tail -n 1)
echo "dummyfilename found to be "$dummyfilename
runnames=$(head -n 4 $masterfile | tail -n 1)
echo "Runnames found to be "$runnames
var_to_change=$(head -n 3 $masterfile | tail -n 1)
echo "var_to_change found to be "$var_to_change
changed_vars=$(head -n 5 $masterfile | tail -n 1)
echo "changed_vars found to be "$changed_vars
current_dir=$(pwd)
echo "current_dir is found to be "$current_dir
echo ''

# Loop through, replacing lines in the dummy_paramfile to create N new paramfiles :)
no_runs=$(echo $runnames | tr ',' '\n' | wc -l)
echo "no_runs found to be "$no_runs

for i in $(seq 1 $no_runs);
do
i=$(echo $i | tr -d '[:space:]')
runname=$(echo $runnames | tr ',' '\n' | head -n $i | tail -n 1)
changed_var=$(echo $changed_vars | tr '|' '\n' | head -n $i | tail -n 1)
echo '	Forming paramfile for run'$runname
cp $dummyfilename $runname.par
sed -i .bak -e "s/run_name=.*/run_name=$runname/" $runname.par
sed -i .bak -e "s/$var_to_change=.*/$var_to_change=$changed_var/" $runname.par
sed -i .bak -e "s:output_folder=.*:output_folder=$current_dir:" $runname.par
rm *.bak
done
echo "Done."

# Now run the first paramfile.......
if test -f "/Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/execute_model.py"; 
then 
echo "WE ARE WORKING IN MAC MODE"
for i in $(echo $runnames | tr ',' '\n')
do
echo $i
cmd="python /Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/execute_model.py $i.par"
echo $cmd
xterm -hold -e "source /Users/mlees/anaconda3/etc/profile.d/conda.sh; conda activate pygmt; python /Users/mlees/Documents/RESEARCH/Land_Subsidence/Local_Scale/MODEL/execute_model.py $i.par" &
read -t 5 -p "I am going to wait for 5 seconds only ..."
done
fi


if test -f "/home/mlees/Land_Subsidence/Local_Scale/MODEL/execute_model.py";
then
echo "WE ARE WORKING IN Linux MODE"
for i in $(echo $runnames | tr ',' '\n')
do
echo $i
cmd="python /home/mlees/Land_Subsidence/Local_Scale/compaction_model/execute_model.py $i.par"
echo $cmd
xterm -hold -e "source /home/mlees/anaconda3/etc/profile.d/conda.sh; conda activate compaction_model; python /home/mlees/Land_Subsidence/Local_Scale/compaction_model/execute_model.py $i.par" &
disown
read -t 50 -p "I am going to wait for 50 seconds now ..."
done
fi
