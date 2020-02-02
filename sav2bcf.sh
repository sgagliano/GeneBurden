fileprefix=$1
chrom=$2

##convert SAV format to BCF
##only extract variants of interest for your burden test:
##`--targets-file chr${chrom}.LoF.regions.txt` to speed up conversion (see end of script for example input)

/net/wonderland/home/lefaivej/savvy/build/sav export -d DS ${fileprefix}.sav | bcftools view --targets-file chr${chrom}.LoF.regions.txt -O b -o ${fileprefix}.bcf

tabix -f ${fileprefix}.bcf

##chr$1.LoF.regions.txt
#chr1    69602   69802
#chr1    930123  944135
#chr1    944993  959138
