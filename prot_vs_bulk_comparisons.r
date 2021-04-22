install.packages("ggpubr")

#get bulk rnaseq pseudotimes and see if proteomic-derived pseudotimes are correlated with rna-seq pseudotimes, for those patients
#included in both studies (n=107 in females)

#firt run protein_lineage_monocle2.R to get Fvariables file, or read in variables file from data folder if saved
read.csv(file="~/prot-lineage/data/Fvariables")
#also run rerun_rnaseq_lineage.r or use saved gene list from that script:
read.csv(file="~/prot-lineage/data/genes_rnaseq")

#upload pseudotimes from bulk rna-seq lineage paper
p <- synapser::synGet('syn22822695')
rnaseq_pstimes <- read.csv(p$path)
rnaseq_pstimes$rnaseq_pseudotime_sc <- scale(rnaseq_pstimes$Pseudotime, center=F)
names(rnaseq_pstimes)[names(rnaseq_pstimes) == "SampleID"] <- "rnaseqID"
names(rnaseq_pstimes)[names(rnaseq_pstimes) == "State"] <- "rnaseq_State"
names(rnaseq_pstimes)[names(rnaseq_pstimes) == "Pseudotime"] <- "rnaseq_Pseudotime"

corrs <- merge(rnaseq_pstimes, Fvariables, by="rnaseqID")
cor(corrs$rnaseq_pseudotime_sc, corrs$pseudotime_sc, method="pearson")
cor(corrs$rnaseq_pseudotime_sc, corrs$pseudotime_sc, method="spearman")
plot(corrs$rnaseq_pseudotime_sc, corrs$pseudotime_sc)
abline(lm(corrs$rnaseq_pseudotime_sc~corrs$pseudotime_sc))
lines(lowess(corrs$rnaseq_pseudotime_sc,corrs$pseudotime_sc))

tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="Pseudotime", y="rnaseq_Pseudotime", 
                  add="reg.line",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_diagnosis.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="Pseudotime", y="rnaseq_Pseudotime", color = "diagnosis",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_braak.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="Pseudotime", y="rnaseq_Pseudotime", color = "braaksc",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_cerad.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="Pseudotime", y="rnaseq_Pseudotime", color = "ceradsc",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_cogdx.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="Pseudotime", y="rnaseq_Pseudotime", color = "cogdx",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()









