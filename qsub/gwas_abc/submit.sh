# args1: middle name of config file
# args2 (optional): preprocessing script (to generate GWAS table)

mkdir -p logs

if [[ ! -z $2 ]]
then
  qsub -v CONFIG=$1,PREP=$2 run.qsub
else
  qsub -v CONFIG=$1 run.qsub
fi
