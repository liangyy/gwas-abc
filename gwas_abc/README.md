This is the main functionality of this repo.
Given the leading variants of the GWAS, we do the following to obtain GWAS signal to gene map:

1. Find ABC enhancers that overlaps with the LD block of the leading variants.
2. Grab the genes that these ABC enhancers are linked to. 
3. Summarize the results for each leading variant by listing the genes along with the ABC score.

Input files:

* GWAS leading variants.
* LD block definition.
* A list of ABC samples.


Test command
```
python gwas_abc.py \
  --gwas_leading_variants test_gwas.tsv.gz chromosome position snpid trait \
  --ld_block /home/t.cri.yliang/labshare/mv_from_scratch/repo_new/rotation-at-imlab/analysis/allelic_heterogeneity/data/ld_independent_regions.txt \
  --abc_pattern /scratch/t.cri.yliang/ABC-GWAS/ABC-by-biosample/AllPredictions.AvgHiC.ABC0.02.ModelRegions.{biosample}.cleanup.tsv.gz 7 21 \
  --abc_sample_yaml test_abc.yaml \
  --bedtools_pylib ~/labshare/GitHub/misc-tools/bedtools_pylib \
  --liftover /gpfs/data/im-lab/nas40t2/yanyul/data/hg38ToHg19.over.chain.gz /gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/liftover_snp \
  --output test.gz
```
