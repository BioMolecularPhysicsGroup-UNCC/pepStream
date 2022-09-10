#!/bin/bash

module load anaconda3

cd ./Filters
tail --lines 2 ./*/output* > temp_output.txt
rm topTenPercent_cMoRFs.txt
rm theBottomTenPercent_cMoRFs.txt
cat contigous_seq_cMoRFs_*.fasta > contigous_seq_cMoRFs.fasta
cp ../format_output.py .
# might need to create the full fasta list here
python3 format_output.py --f1 temp_output.txt --cmorfs contigous_seq_cMoRFs.fasta

# dir with zDockSpecific models > zDockModels/specific_models

while read file; do
	echo "$file in file..."
        tempName=$file  
        cd $tempName*
	cp -r ../zDockModels .
        cp ../run_specificity.sh .
	chmod 775 *
        sbatch run_specificity.sh
        cd ..
done < topTenPercent_cMoRFs.txt
#done < theBottomTenPercent_cMoRFs.txt

# run after running this to work up specific study
# tail --lines 4 ./*/specific_output* > specific_temp_output.txt
# rm specific_cMoRFs_output.txt
# python3 format_output_specific.py --f1 specific_temp_output.txt
