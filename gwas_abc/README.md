This is the main functionality of this repo.
Given the leading variants of the GWAS, we do the following to obtain GWAS signal to gene map:

1. Find ABC enhancers that overlaps with the LD block of the leading variants.
2. Grab the genes that these ABC enhancers are linked to. 
3. Summarize the results for each leading variant by listing the genes along with the ABC score.

Input files:

* GWAS leading variants.
* LD block definition.
* A list of ABC samples.
