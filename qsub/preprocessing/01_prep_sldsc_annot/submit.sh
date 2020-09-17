#ARGS1: list of biosample
if [[ -z $1 ]]
then
  biosamplelist=/gpfs/data/im-lab/nas40t2/yanyul/ABC_enhancer2gene/AllPredictions.AvgHiC.ABC0.02.ModelRegions.txt.gz.all_celltype_list
else
  biosamplelist=$1
fi

for i in `cat $biosamplelist`
do
  qsub -v BIOSAMPLE=$i run.qsub
done
