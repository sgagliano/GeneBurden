# GeneBurden
From minimac4 BCF/VCF.gz input create a gene burden VCF.gz file

The gene burden is the sum of alternate allele dosages for a specified set of variants (e.g. LoF variants within a gene)

### Example: 

Create a gene burden VCF file from variants in the BED file (2 columns expected: 1st variant ID that matches variant ID column in inGZVCF, 2nd: gene/region name)
Include rare variants with an imputation R2 of at least 0.1 (`--minR2 0.1`) and an alternate allele frequency of less than 0.005 (`--maxAF 0.005`)

`python AltAlleleBurden-vcf-SUMMATION.py --inGZVCF ${vcffile} --outGZAnno chr${chrom}.Burden.vcf.gz --minR2 0.1 --maxAF 0.005 --BED chr${chrom}.LoF.list.txt`
