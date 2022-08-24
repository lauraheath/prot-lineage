#look at correlation between LFQ log2FC and TMT log2FC
library(ggpubr)
#compare LFQ agora data and TMT agora data results

p <- synapser::synGet('syn18689335')
lfq_mat <- read.csv(p$path)

tmt_mat <- allresults

lfq_dlpfc <- subset(lfq_mat, lfq_mat$Tissue=='DLPFC', select=c(UniqID, Log2_FC))
lfq_ant <- subset(lfq_mat, lfq_mat$Tissue=='AntPFC', select=c(UniqID, Log2_FC))
lfq_mfg <- subset(lfq_mat, lfq_mat$Tissue=='MFG', select=c(UniqID, Log2_FC))
lfq_tcx <- subset(lfq_mat, lfq_mat$Tissue=='TCX', select=c(UniqID, Log2_FC))

names(lfq_dlpfc)[names(lfq_dlpfc) == "Log2_FC"] <- "dlpfc_Log2FC"
names(lfq_ant)[names(lfq_ant) == "Log2_FC"] <- "ant_Log2FC"
names(lfq_mfg)[names(lfq_mfg) == "Log2_FC"] <- "mfg_Log2FC"
names(lfq_tcx)[names(lfq_tcx) == "Log2_FC"] <- "tcx_Log2FC"
names(tmt_mat)[names(tmt_mat) == "Log2_FC"] <- "tmt_Log2FC"

tmt_mat <- dplyr::left_join(tmt_mat, lfq_dlpfc)
tmt_mat <- dplyr::left_join(tmt_mat, lfq_ant)      
tmt_mat <- dplyr::left_join(tmt_mat, lfq_mfg)        
tmt_mat <- dplyr::left_join(tmt_mat, lfq_tcx)   


cor(tmt_mat$dlpfc_Log2FC,tmt_mat$tmt_Log2FC, method="pearson", use="complete.obs")


ggscatter(tmt_mat, x = "dlpfc_Log2FC", y = "tmt_Log2FC", use="complete.obs",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lfq dlpfc", ylab = "TMT")

ggscatter(tmt_mat, x = "dlpfc_Log2FC", y = "ant_Log2FC", use="complete.obs",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lfq dlpfc", ylab = "lfq antPFC")

ggscatter(tmt_mat, x = "dlpfc_Log2FC", y = "mfg_Log2FC", use="complete.obs",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lfq dlpfc", ylab = "lfq mfg")

ggscatter(tmt_mat, x = "dlpfc_Log2FC", y = "tcx_Log2FC", use="complete.obs",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lfq dlpfc", ylab = "lfq tcx")


ggscatter(tmt_mat, x = "ant_Log2FC", y = "tmt_Log2FC", use="complete.obs",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lfq ant", ylab = "TMT")
ggscatter(tmt_mat, x = "mfg_Log2FC", y = "tmt_Log2FC", use="complete.obs",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lfq mfg", ylab = "TMT")
ggscatter(tmt_mat, x = "tcx_Log2FC", y = "tmt_Log2FC", use="complete.obs",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "lfq tcx", ylab = "TMT")
