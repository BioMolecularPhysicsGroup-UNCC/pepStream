peptide=M-P10
echo "Summary of model Zdock binding" > ${peptide}_summary.out
for output in *.out
do

echo ${output} >> ${peptide}_summary.out
sed -n '6,10p' ${output} >> ${peptide}_summary.out

done
