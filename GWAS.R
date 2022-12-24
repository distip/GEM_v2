
library(rMVP)
# Prepare data for rMVP, will convert files in a big matrix format. limitation will be the size storage and not RAM
MVP.Data(fileVCF ="BGEM.all.imputed.vcf", # Genotype in hapmap format. inds: 329        markers:14.7 M
         filePhe="pheno.csv", # Phenotype, the first column is taxa name, the subsequent columns are traits
         sep.hmp="\t",
         sep.phe=",",
         SNP.effect="Add",
         fileKin=TRUE,
         filePC=TRUE,
         #priority="memory", #"speed" or "memory"
         #maxLine=10000, #number of SNPs, only used for saving memory when calculate kinship matrix
         out="/work/schnablelab/deniz/BGEM/GWAS/analysis/mvp.hmp") # prefix of output file name, will be placed in the working directory
genotype = attach.big.matrix("mvp.hmp.geno.desc")
phenotype = read.table("mvp.hmp.phe",head=TRUE)
map =  read.table("mvp.hmp.geno.map" , head = TRUE)
imMVP <- MVP(
  phe=phenotype[,c(1,550)], #1 = genotype, 8 = trait
  geno=genotype,
  map=map,
  #K=Kinship,
  #CV.GLM=Covariates,     ##if you have additional covariates, please keep there open.
  #CV.MLM=Covariates,
  #CV.FarmCPU=Covariates,
  #nPC.GLM=5,      ##if you have added PC into covariates, please keep there closed.
  nPC.MLM=3,
  nPC.FarmCPU=3,
  priority="speed",       ##for Kinship construction
  #ncpus=10,
  vc.method="EMMA",      ##only works for MLM
  maxLoop=10,
  method.bin="FaST-LMM",      ## "EMMA", "FaST-LMM", "static" (#only works for FarmCPU)
  #permutation.threshold=TRUE,
  #permutation.rep=100,
  threshold=0.05,
  method = c("MLM","FarmCPU") #method=c("GLM", "MLM", "FarmCPU")