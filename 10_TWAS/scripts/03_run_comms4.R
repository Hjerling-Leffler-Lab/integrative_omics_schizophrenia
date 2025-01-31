# run CoMM-S4 (pmid=34616426) for brain cell type eQTLs and scz2022-EUR

library(CoMM)

args <- commandArgs(trailingOnly = TRUE)
i=as.numeric(args)[1]
CT=as.character(args)[2]
lvl=as.character(args)[3]

print(i)
print(CT)
print(lvl)

INP="./10_TWAS/data/"
REF=paste("./workflow/REF/perCHR_",lvl,"_eqtlp0.1/",sep="")
OUTDIR="./10_TWAS/data/TWAS_output/"

stringname1 = paste(INP,lvl,"/",CT,"/INPUT1_chr",i,".eQTL_combined_p0.1.txt",sep="") 
stringname2 = paste(INP,lvl,"/",CT,"/INPUT2_chr",i,"_gwas_scz2022.txt",sep="") 
stringname3 = paste(REF,"HRC.r1-1.EGA.GRCh37.noATCG.nomhc_eqtlp0.1_chr",i,sep="") 
stringname4 = paste(INP,lvl,"/",CT,"/INPUT4_chr",i,".eQTL_combined_p0.1",sep="")
stringname5 = paste(REF,"HRC.r1-1.EGA.GRCh37.noATCG.nomhc_eqtlp0.1_chr",i,sep="")
px = 1
lam = 0.95;
coreNum = as.numeric(args)[4];

print(stringname1)
print(stringname2)
print(stringname3)
print(stringname4)
print(stringname5)

fm=CoMM_S4_testing_mt(stringname1,stringname2,stringname3,stringname4,stringname5,px,lam,coreNum)

OUTPUT=paste(OUTDIR,lvl,"/",CT,"/OUTPUT_chr",i,"_",CT,"_scz2022.RData",sep="")
save.image(file=OUTPUT)
