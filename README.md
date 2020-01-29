# GeneBurden
From minimac4 BCF/VCF.gz input, create a region/gene burden VCF.gz file

The region/gene burden in this script is defined as the sum of alternate allele dosages for a specified set of variants (e.g. putative loss of function variants, LoF within each region/gene)

### Example: 

##### Create a region/gene burden VCF.gz file

Create a gene burden VCF.gz file from variants in the BED file (2 columns expected: 1st variant ID that matches variant ID column in inGZVCF, 2nd: gene/region name)
Include rare variants with an imputation R2 of at least 0.1 (`--minR2 0.1`) and an alternate allele frequency of less than 0.005 (`--maxAF 0.005`)

`python AltAlleleBurden-vcf-SUMMATION.py --inGZVCF ${vcffile} --outGZAnno chr${chrom}.Burden.vcf.gz --minR2 0.1 --maxAF 0.005 --BED chr${chrom}.LoF.list.txt`

Hint: to create a `--BED` file of LoF variants (as defined by VEP) refer to `extract_lof_from_vep.py`

##### Test the regions/genes for association with your trait of interest

See the `SAIGE/` directory for scripts for ```step1``` (compute the GRM and variance ratio using a set of high-quality pruned SNPs) and ```step2``` (input the output from `--outGZAnno` to test each region/gene for association with the phenotype of interest)

For more information refer to SAIGE on GitHub: https://github.com/weizhouUMICH/SAIGE and the publication Zhou _et al._ Nature Genetics 50, 1335–1341 (2018)

