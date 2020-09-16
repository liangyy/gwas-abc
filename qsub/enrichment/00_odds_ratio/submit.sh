#ARGS1: list of GWAS
gwaslist=../../../trait_list.txt
nperbatch=10
if [[ ! -d gwas_list ]]
then
  mkdir -p gwas_list
  split -l $nperbatch -d $gwaslist gwas_list/batch
fi

for i in `ls gwas_list/*`
do
  j=`echo $i | sed 's#gwas_list/##g'`
  qsub -v BATCH=$j run.qsub
done
