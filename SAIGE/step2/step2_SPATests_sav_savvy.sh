#Example run:
#sh step2_SPATests_sav_savvy.sh  X174 chr1 1 10822 19868066
#where X174 is our trait name, chr1 is the chunk name, and we test variants within 1:10822:19868066

trait=$1 #X..
chunkName=$2 #your burden file name before the .Burden.vcf.gz ending
chrom=$3 #1,2,...22,X
start=$4
end=$5

Rscript step2_SPATests.R \
        --vcfFile= ${chunkName}.Burden.vcf.gz \
        --vcfFileIndex=${chunkName}.Burden.vcf.gz.tbi \
        --chrom=chr${chrom} \
        --start=$start \
        --end=$end \
        --GMMATmodelFile=../step1/output/${trait}.rda \
        --varianceRatioFile=../step1/output/${trait}.varianceRatio.txt \
        --numLinesOutput=1000 \
        --sampleFile=vcfsample.list \
        --SAIGEOutputFile=output/${trait}.${chunkName}.SAIGE.txt \
        --IsOutputAFinCaseCtrl=TRUE

