module load gcc/6.2.0
module load R/3.6.3

filename=OMIM-LD-block.tsv.gz
if [[ -f $filename ]]
then
  wget \
    -O $filename.tmp \
    https://bitbucket.org/yanyul/rotation-at-imlab/raw/a4ad9182aac280b1e4b84c62b4473acdcb56866a/analysis/fdr_power_specificity/companion_figure_finalized/summary_on_expression_cleanup/logistic-based-test.datamatrix-with-gene-info.OMIM-LD-block-PrediXcan-MASH-EUR.tsv

  Rscript -e "f = read.table('OMIM-LD-block.tsv.gz.tmp', sep = '\t', header = T); \
  f = f[ f$is_omim];\
  chr = unlist(lapply(strsplit(f$lead_var), function(x) {x[1]}));\
  pos = unlist(lapply(strsplit(f$lead_var), function(x) {x[2]}));\
  out = data.frame(chromsome = chr, position = pos, lead_var = f$lead_var, trait = f$trait);\
  gz = gzfile('OMIM-LD-block.tsv.gz', 'w');\
  write.table(gz, sep='\t', row=F, col=T, quo=F);\
  close(gz)
  "
  rm $filename.tmp
fi