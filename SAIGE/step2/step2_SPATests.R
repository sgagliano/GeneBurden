options(stringsAsFactors=F)

library(SAIGE, lib.loc="/net/snowwhite/home/sarahgag/R/x86_64-pc-linux-gnu-library/3.3") ##INSERT YOUR PATH TO SAIGE HERE
library(optparse)


option_list <- list(
  make_option("--savFile", type="character",default="",
    help="path to savFile."),
  make_option("--savFileIndex", type="character",default="",
    help="path to savFile index file."),
  make_option("--chrom", type="character",default="",
    help="chromosome."),
  make_option("--start", type="numeric",default=0,
    help="start position."),
  make_option("--end", type="numeric",default=200000000000,
    help="end position."),
  make_option("--sampleFile", type="character",default="",
    help="File contains one column for IDs of samples in the bgen file, no header"),
  make_option("--GMMATmodelFile", type="character",default="",
    help="path to the input file containing the glmm model"),
  make_option("--varianceRatioFile", type="character",default="",
    help="path to the input file containing the variance ratio"),
  make_option("--SAIGEOutputFile", type="character", default="",
    help="path to the output file containing the SAIGE test results"),
  make_option("--numLinesOutput", type="numeric", default=10000,
    help="output every nth markers"),
  make_option("--IsOutputAFinCaseCtrl", type="logical", default=TRUE,
    help="outputAFinCaseCtrl [default='TRUE']")
)


parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

try(if(length(which(opt == "")) > 0) stop("Missing arguments"))


SPAGMMATtest(savFile=opt$savFile,
             savFileIndex=opt$savFileIndex,
             chrom=opt$chrom,
             start=opt$start,
             end=opt$end,       
             sampleFile=opt$sampleFile,
             GMMATmodelFile=opt$GMMATmodelFile,
             varianceRatioFile=opt$varianceRatioFile,
             SAIGEOutputFile=opt$SAIGEOutputFile,
             numLinesOutput=opt$numLinesOutput,
             IsOutputAFinCaseCtrl=TRUE


