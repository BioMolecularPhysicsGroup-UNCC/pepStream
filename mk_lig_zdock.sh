#!/bin/sh
mkdir output_model1_151119
mkdir output_model2_151119
mkdir output_model3_151119
mkdir output_model4_151119
mkdir output_model5_151119

for VARIABLE in 1 2 3 4 5 6 7 8 9 10
do
create.pl model1_DELLA_Model${VARIABLE}.out 5
for nam in complex.*
do
    newname=${nam%.pdb}model1_DELLA${VARIABLE}.pdb
    mv $nam $newname
done
mv complex* ./output_model1_151119/

create.pl model2_DELLA_Model${VARIABLE}.out 5
for nam in complex.*
do
    newname=${nam%.pdb}model2_DELLA${VARIABLE}.pdb
    mv $nam $newname
done
mv complex* ./output_model2_151119/

create.pl model3_DELLA_Model${VARIABLE}.out 5
for nam in complex.*
do
    newname=${nam%.pdb}model3_DELLA${VARIABLE}.pdb
    mv $nam $newname
done
mv complex* ./output_model3_151119/

create.pl model4_DELLA_Model${VARIABLE}.out 5
for nam in complex.*
do
    newname=${nam%.pdb}model4_DELLA${VARIABLE}.pdb
    mv $nam $newname
done
mv complex* ./output_model4_151119/

create.pl model5_DELLA_Model${VARIABLE}.out 5
for nam in complex.*
do
    newname=${nam%.pdb}model5_DELLA${VARIABLE}.pdb
    mv $nam $newname
done
mv complex* ./output_model5_151119/
done
