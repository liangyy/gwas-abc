#ARGS1: list of biosample
if [[ -z $1 ]]
then
  gwaslist=../../../trait_list.txt
else
  gwaslist=$1
fi

for i in `cat $gwaslist`
do
  qsub -v GWAS=$i run.qsub
done
