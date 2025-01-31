# readme.txt
# for folder: 9_sc_eQTL/
# by Shuyang Yao (shuyang.yao@ki.se) 2024-11-06
# Goal: run eQTL analysis at cell type level; 3CT: 3 cell classes, namely Excitatory, Inhibitory, and NonNeuronal; 16CT: 16 cell clusters, see paper.
# folder structure below:
|__ scripts/:
01_a_prepare_vcf.Rmd:        # prepare genotype vcf input to QTLtools. The files are stored in EGA.
01_b_prepare_expression_input.Rmd:        # prepare phenotype input to QTLtools.
01_c_prepare_covariates.Rmd:        # prepare covariate input to QTLtools.
02_eQTL_16CT.Rmd:        # run eQTL analysis in QTLtools for 16 cell type level.
02_eQTL_3CT.Rmd:        # run eQTL analysis in QTLtools for 3 cell type level.
rscript_01_b1_Filter_Genes_16CTs.R:        # R script called in 01_b_prepare_expression_input.Rmd, filter genes at 16 cell type level.
rscript_01_b1_Filter_Genes_3CTs.R:        # R script called in 01_b_prepare_expression_input.Rmd, filter genes at 3 cell type level.
rscript_01_b2_Get_Exp_Input_16CTs.R:        # R script called in 01_b_prepare_expression_input.Rmd, get expression input file at 16 cell type level.
rscript_01_b2_Get_Exp_Input_3CTs.R:        # R script called in 01_b_prepare_expression_input.Rmd, get expression input file at 3 cell type level.
rscript_01_c_ExpPC.R:        # R script called in 01_c_prepare_covariates.Rmd
rscript_02_eQTL_qvalue_16CT.R:        # R script called in 02_eQTL_16CT.Rmd, get qvalue for eQTL output
rscript_02_eQTL_qvalue_3CT.R:        # R script called in 02_eQTL_3CT.Rmd, get qvalue for eQTL output
|__ output/        # eQTL results, filtered for q-value≤0.05
       |__ 3CT/: ${CellType}_all.chr.eQTL_combined_qvalue0.05.tsv.gz # eQTL sumstats at cell class level. Will share when paper is published
       |__ 16CT/: ${CellType}_all.chr.eQTL_combined_qvalue0.05.tsv.gz # eQTL sumstats at cell type level. Will share when paper is published
       |__ plot folder for each relevant plot panel

#-- end --#

