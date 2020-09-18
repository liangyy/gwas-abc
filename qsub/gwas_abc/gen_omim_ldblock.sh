module load gcc/6.2.0
module load R/3.6.3

filename=OMIM-LD-block.tsv.gz
if [[ ! -f $filename ]]
then
  wget \
    -O $filename.tmp \
    https://bitbucket.org/yanyul/rotation-at-imlab/raw/a4ad9182aac280b1e4b84c62b4473acdcb56866a/analysis/fdr_power_specificity/companion_figure_finalized/summary_on_expression_cleanup/logistic-based-test.datamatrix-with-gene-info.OMIM-LD-block-PrediXcan-MASH-EUR.tsv

  Rscript gen_omim_ldblock.R 
  rm $filename.tmp
fi
