#ARGS1: run script name (default run.qsub)

gwaslist=../../../trait_list.txt
nperbatch=10
if [[ ! -d gwas_list ]]
then
  mkdir -p gwas_list
  split -l $nperbatch -d $gwaslist gwas_list/batch
fi

if [[ -z $1 ]]
then 
  runscript=run.qsub
else
  runscript=$1
fi

for i in `ls gwas_list/*`
do
  j=`echo $i | sed 's#gwas_list/##g'`
  qsub -v BATCH=$j $runscript
done
