#ARGS1: list of biosample
gwaslist=../../../trait_list.txt

for i in `cat $gwaslist`
do
  qsub -v GWAS=$i run.qsub
done