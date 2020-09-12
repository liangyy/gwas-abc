# gwas-abc

A ready-to-go pipeline to utilize ABC models proposed [here](https://www.nature.com/articles/s41588-019-0538-0) and calculated [here](https://www.biorxiv.org/content/10.1101/2020.09.01.278093v1) to GWAS summary statistics.

The pipeline do two major things:

* Perform biosample enrichment and S-LDSC regression to examine the enrichment of enhancer for each biosample.
* Link GWAS loci to gene via ABC enhancer-gene map.
