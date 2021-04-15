##example run:
#Rscript step1_spagmmat-BINARY.R X174 
#Where X174 is the name of the trait in our ydatFile that we want to test for association

#this script performs step1 of SAIGE (computing GRM and variance ratio)
#the input for step1 (<mainPath>/step1/output/trait.rda, <mainPath>/step1/output/varianceRatio.txt) will be used in step2

rm(list=ls())

options(echo=TRUE, stringsAsFactors=F) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
trait = as.character(args[1])
 
cat("trait", trait)

library(Rcpp)
library(SAIGE, lib.loc="/net/snowwhite/home/sarahgag/R/x86_64-pc-linux-gnu-library/4.0") #REPLACE CONTENT IN QUOTES YOUR PATH TO SAIGE
#how I downloaded SAIGE:
#library('devtools')
#devtools::install_github("leeshawn/MetaSKAT") 
#devtools::install_github("weizhouUMICH/SAIGE") 

mainPath = "/net/inpsyght/disk2/sarahgag/singlevariant" #REPLACE CONTENT IN QUOTES TO YOUR MAIN PATH (where step1 and step2 are sub-directories)

ydatFile = paste0("Phenotype_with_covariate_file.txt")

outputFolder=paste0(mainPath, "/step1/output/")

modelOut=paste0(outputFolder,trait)


genoFile = "<INSERT PATH TO PLINK FILES, OMIT THE BED/BIM/FAM SUFFIX WITH PRUNED SET OF GENOME-WIDE HIGH QUALITY VARIANTS>"


fitNULLGLMM(plinkFile = genoFile,
                phenoFile = ydatFile,
                phenoCol = trait,
                traitType = "binary",
                invNormalize = FALSE,
                covarColList = c("Sex", "birthYear","PC1.wb","PC2.wb","PC3.wb","PC4.wb"), #INSERT YOUR COVARIATE NAMES HERE
                qCovarCol = NULL,
                sampleIDColinphenoFile = "IID", #INSERT COLUMN NAME WITH SAMPLE IDS in ydatFile
                tol=0.02,
                maxiter=20,
                tolPCG=1e-5,
                maxiterPCG=500,
                nThreads = 32,
                #Cutoff = 2,
                numMarkers = 100,
                skipModelFitting = FALSE,
                outputPrefix = modelOut,
		memoryChunk = 2)
