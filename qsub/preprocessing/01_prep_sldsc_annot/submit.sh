#ARGS1: list of biosample
biosamplelist=/gpfs/data/im-lab/nas40t2/yanyul/ABC_enhancer2gene/AllPredictions.AvgHiC.ABC0.02.ModelRegions.txt.gz.all_celltype_list

for i in `cat $biosamplelist`
do
  qsub -v BIOSAMPLE=$i run.qsub
done