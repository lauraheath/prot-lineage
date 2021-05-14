Log2_Normalized <- readRDS(file="~/prot-lineage/data/Log2_Normalized.rds")
Meta <- readRDS(file="~/prot-lineage/data/Meta.rds")

#function to run monocle analysis
RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                                  gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  #HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = ~pmi+educ)
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  #HSMM <- orderCells(HSMM, reverse=TRUE)
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}

#run get_proteins_script_JG.r to get metadata and log2 normalized protein matrix
synapser::synLogin()

Dat <- Log2_Normalized
Dat[is.na(Dat)] <- 0

#get list of proteins that are differentially expressed between AD case & control (Jake's analysis)
p <- synapser::synGet('syn24216817')
ADgenes_prot <- read.csv(p$path)
ADgenes_prot1 <- subset(ADgenes_prot, ADgenes_prot$PVal<0.05)
dim(ADgenes_prot1)



genes2<-c()
for (gene in unique(c(as.vector(ADgenes_prot1$Peptide)))){
  if (gene %in% rownames(Dat)){
    genes2 <- c(genes2,which(rownames(Dat)==gene))
  }
}
length(genes2)
Dat2 <- Dat[genes2,]
dim(Dat2)



#Keeping only female data (msex==0 is female, msex==1 is male; run separately for sex-specific analysis)
In_S <- which(Meta$msex == 0)
#In_S <- which(Meta$msex == 1)
Dat2 <- Dat2[,In_S]
Meta2 <- Meta[In_S,]

gene_short_name <- rownames(Dat2)
temp <- Dat2
temp2 <- Meta2


temp2$braaksc <- factor(temp2$braaksc,levels = c(0:6))
temp2$ceradsc <- factor(temp2$ceradsc,levels = c(1:4))
temp2$cogdx <- factor(temp2$cogdx,levels = c(1:6))

rownames(temp)<-NULL
rownames(temp2)<-NULL

#gene_df <- as.data.frame(gene_short_name)
#write.csv(gene_df, file="~/prot-lineage/data/prot_genes.csv", row.names=FALSE)

MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)

saveRDS(MonRun, file="~/prot-lineage/data/FEMALE_monocleobject_matrix.RDS")
#saveRDS(MonRun, file="~/prot-lineage/data/MALE_monocleobject_matrix.RDS")


g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g


g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g<- plot_cell_trajectory(MonRun,color_by = "batch",show_branch_points=F,use_color_gradient = F,cell_size = 1.5) + theme(legend.position="none")
g <- g + ggplot2::scale_color_viridis_d() + theme(legend.position="none")
g <- g + ggplot2::labs(color="diagnosis")
g

#collapse some states due to small numbers and rearrange
table(MonRun$State)

#female samples
# MonRun$State2 <- MonRun$State
# MonRun$State2[MonRun$State == 3] <- 2
# MonRun$State2[MonRun$State == 4] <- 3
# MonRun$State2[MonRun$State == 5] <- 4
# MonRun$State2[MonRun$State == 6] <- 5
# MonRun$State2[MonRun$State == 7] <- 6

#male samples:
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 3] <- 2
MonRun$State2[MonRun$State == 11] <- 3
MonRun$State2[MonRun$State == 7] <- 6
MonRun$State2[MonRun$State == 8] <- 7
MonRun$State2[MonRun$State == 10] <- 7

MonRun$State2 <- as.numeric(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g


#saveRDS(MonRun, file="~/prot-lineage/data/Female_monocleObject.rds")
saveRDS(MonRun, file="~/prot-lineage/data/Male_monocleObject.rds")

###### FIGURES #########
#tiff(file='~/prot-lineage/figures/FEMALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "braaksc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Braak Score")
g
dev.off()

MonRun$APOE <- factor(MonRun$APOE,levels=c(0,1,2))
#tiff(file='~/prot-lineage/figures/FEMALE_tree_apoe.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_tree_apoe.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "APOE",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="APOE e4 Dosage")
g
dev.off()


#tiff(file='~/prot-lineage/figures/FEMALE_tree_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_tree_cerad.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "ceradsc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="CERAD Score")
g
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_tree_cogdx.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_tree_cogdx.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "cogdx",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Cognitive Diagnosis")
g
dev.off()


#tiff(file='~/prot-lineage/figures/FEMALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
dev.off()


#tiff(file='~/prot-lineage/figures/FEMALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=scale(Pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=cogdx, y=scale(Pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Cognitive\nDiagnosis",y="Pseudotime",x="Cognitive Diagnosis")
g
dev.off()





######## for stats & other figures, create a dataframe with all relevant covariates & pseudotime & state

x <- list()
x$SampleID <- MonRun$batchChannel
x$rnaseqID <- MonRun$rnaseq_id
x$State2 <- MonRun$State2
x$Pseudotime <- MonRun$Pseudotime
x$diagnosis <- MonRun$diagnosis
x$braaksc <- MonRun$braaksc
x$ceradsc <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx
x$apoe <- MonRun$APO
x$educ   <- MonRun$educ
x$pmi <- MonRun$pmi
x$batch <- MonRun$batch
x$mmse <- MonRun$cts_mmse30_lv
x$age_death <- MonRun$age_death
x$rna_seq_sample <- MonRun$rnaseq
x$SampleID <- as.character(x$SampleID)
Fvariables <- as.data.frame(x)
Fvariables$pseudotime_sc <- scale(Fvariables$Pseudotime, center=F)

casecontrolF <- subset(Fvariables, Fvariables$diagnosis=='AD'|Fvariables$diagnosis=='Control')
casecontrolF$diag2 <- ifelse(casecontrolF$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ pseudotime_sc,casecontrolF,family='binomial'))
#tiff(file='~/prot-lineage/figures/FEMALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)

g <- ggplot(casecontrolF,aes(x=diagnosis,
                        y=pseudotime_sc,
                        color=diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g
dev.off()

braakfit <- MASS::polr(braaksc ~ pseudotime_sc,Fvariables)
ceradfit <- MASS::polr(ceradsc ~ pseudotime_sc,Fvariables)
cogdxfit <- MASS::polr(cogdx ~ pseudotime_sc,Fvariables)

cat('braak p-value: ',pt(abs(summary(braakfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
cat('cerad p-value: ',pt(abs(summary(ceradfit)$coef[1,3]),ceradfit$df.residual,lower.tail=F)*2,'\n')
cat('cogdx p-value: ',pt(abs(summary(cogdxfit)$coef[1,3]),cogdxfit$df.residual,lower.tail=F)*2,'\n')



#look for correlations with GWAS LOAD genes
ad_gwas <- c("CR1",
             "BIN1",
             "INPP5D",
             "HLA-DRB1",
             "TREM2",
             "MEF2C",
             "NME8",
             "CD2AP",
             "NYAP1",
             "EPHA1",
             "PTK2B",
             "CLU",
             "SPI1",
             "MS4A2",
             "PICALM",
             "SORL1",
             "FERMT2",
             "SLC24A4",
             "ABCA7",
             "APOE",
             "CASS4",
             "ECHDC3",
             "ACE",
             "NDUFAF6",
             "ECHDC3",
             "ADAMTS20",
             "SPPL2A",
             "ADAM10",
             "IQCK",
             "MIR142",
             "ACE",
             "ADAMTS1",
             "SUCLG2P4",
             "FST",
             "OARD1",
             "WWOX",
             "MAF",
             "CD55",
             "YOD1",
             "HLA-DRB1",
             "PSMB8",
             "C4A",
             "GPSM3",
             "HLA-DPA1",
             "HLA-DQA1",
             "HLA-DRA",
             "HLA-DRB5",
             "PSMB9",
             "CD2AP",
             "AGFG2",
             "PILRA",
             "EPHB4",
             "C7orf43",
             "GAL3ST4",
             "ZKSCCAN1",
             "FAM131B",
             "PSMC3",
             "ACP2",
             "C1QTNF4",
             "CELF1",
             "MTCH2",
             "NDUFS3",
             "NUP160",
             "MS4A6A",
             "MS4A7",
             "MS4A4A",
             "EED",
             "PICALM",
             "STYX",
             "RIN3",
             "HMHA1",
             "CNN2",
             "WDR18",
             "CASS4")

Dat3 <- Dat2
Dat3$gene_names <- rownames(Dat3)
Dat3$gene_short_name <- gsub("\\|.*", "", Dat3$gene_names)


# dlpfcCPMObj <- synapser::synGet('syn8456638')
# Dat <- data.table::fread(dlpfcCPMObj$path,data.table=F)
sampleIds <- colnames(Dat2)[-120]
sampleIds <- sampleIds[-120]
sampleIds
geneIds <- Dat3$gene_short_name
Dat3$gene_short_name<-NULL
Dat3$gene_names<-NULL
Dat3 <- t(Dat3)
colnames(Dat3) <- geneIds
Dat3 <- data.frame(Dat3,stringsAsFactors=F)
Dat3$sampleId <- sampleIds
dlpfc <- dplyr::left_join(Fvariables,Dat3,by=c('SampleID'='sampleId'))
dlpfc2 <- dlpfc[,17:1566]

corvec <- cor(dlpfc2,dlpfc$Pseudotime,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$geneid %in% ad_gwas,]$cor))

mean(corDfdlpfc$cor)
mean(corDfdlpfc[corDfdlpfc$geneid %in% ad_gwas,]$cor)


corDfdlpfc$adGwas <- corDfdlpfc$geneid %in% ad_gwas
colnames(corDfdlpfc)[3] <- 'LOADGWASGene'
corDfdlpfc$LOADGWASGene2 <- ifelse(corDfdlpfc$LOADGWASGene==FALSE, "NOT GWAS GENE", "GWAS GENE")

#tiff(file='~/prot-lineage/figures/FEMALE_loadgwas_cor.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_loadgwas_cor.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(corDfdlpfc,ggplot2::aes(x=LOADGWASGene2,y=cor,fill=LOADGWASGene2))
g <- g + ggplot2::geom_boxplot() + theme(legend.position="none")
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'LOAD GWAS STATUS',y='Correlation with pseudotime')
g
dev.off()


#save variables file
#write.csv(Fvariables, file="~/prot-lineage/data/Fvariables.csv", row.names=FALSE)
write.csv(Fvariables, file="~/prot-lineage/data/Mvariables.csv", row.names=FALSE)

######### branch-specific differential expression analysis

#need Monrun object, gene_short_name

#there are 6 states
#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))
l2$p_6 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))


for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(MonRun$State2)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$p_3[i] <- tk$s[2,4]
  l2$p_4[i] <- tk$s[3,4]
  l2$p_5[i] <- tk$s[4,4]
  l2$p_6[i] <- tk$s[5,4]
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
  l2$d_6[i] <- tk$s[5,1]
}

#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))

dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)

dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])

df3 <- dplyr::left_join(dfa1,dfb1)

#gene names are currently labeled by gene and uniprot id. need to split:
df4 <- df3
df4$gene_short_name <- gsub("\\|.*", "", df4$gene_names)
#reorder columns
df4 <- df4[,c(1,5,2,3,4)]
#write.csv(df4,file='~/prot-lineage/results/female_anova_stats.csv',quote=F,row.names=F)
write.csv(df4,file='~/prot-lineage/results/male_anova_stats.csv',quote=F,row.names=F)





#run branch-specific enrichment analysis 
library(enrichR)

#use the initial data frame with all p-values & effect sizes from tukey analysis (df2)
dlpfc <- df2
dlpfc$gene_names <- gsub("\\|.*", "", dlpfc$gene_names)
getGeneSets = function(pval,dval,setName,setDescription,summaryMat){
  resPos <- dplyr::filter(summaryMat,.data[[pval]]<0.05, .data[[dval]]>0)
  resNeg <- dplyr::filter(summaryMat,.data[[pval]]<0.05, .data[[dval]]<0)
  res <- list()
  res$pos <- c(paste0(setName,'pos'),paste0(setDescription,'pos'),resPos$gene_names)
  res$neg <- c(paste0(setName,'neg'),paste0(setDescription,'neg'),resNeg$gene_names)
  return(res)
}
runList <- list()
runList$pval <- paste0('p_',2:6)
runList$dval <- paste0('d_',2:6)
runList$setName <- paste0('branch',2:6)
runList$setDescription <- paste0('branchSpecificDEG',2:6)
res <- mapply(getGeneSets,
              runList$pval,
              runList$dval,
              runList$setName,
              runList$setDescription,
              MoreArgs = list(summaryMat=dlpfc),
              SIMPLIFY= FALSE)
res <- unlist(res,recursive = F)

# dummyFun <- function(geneSet){
#   res <- enrichR::enrichr(geneSet[-c(1:2)],databases = c('GO_Biological_Process_2018',
#                                                          'GO_Molecular_Function_2018',
#                                                          'GO_Cellular_Component_2018'))
#   return(res)
# }
# dummyFun2 <- function(listOfGO){
#   listOfGO$GO_Biological_Process_2018 <- dplyr::filter(listOfGO$GO_Biological_Process_2018,Adjusted.P.value <0.05)
#   listOfGO$GO_Molecular_Function_2018 <- dplyr::filter(listOfGO$GO_Molecular_Function_2018,Adjusted.P.value <0.05)
#   listOfGO$GO_Cellular_Component_2018 <- dplyr::filter(listOfGO$GO_Cellular_Component_2018,Adjusted.P.value <0.05)
#   return(listOfGO)
# }

dummyFun <- function(geneSet){
  res <- enrichR::enrichr(geneSet[-c(1:2)],databases = c('Allen_Brain_Atlas_10x_scRNA_2021',
                                                         'Allen_Brain_Atlas_down',
                                                         'Allen_Brain_Atlas_up'))
  return(res)
}
dummyFun2 <- function(listOfGO){
  listOfGO$Allen_Brain_Atlas_10x_scRNA_2021 <- dplyr::filter(listOfGO$Allen_Brain_Atlas_10x_scRNA_2021,Adjusted.P.value <0.05)
  listOfGO$Allen_Brain_Atlas_down <- dplyr::filter(listOfGO$Allen_Brain_Atlas_down,Adjusted.P.value <0.05)
  listOfGO$Allen_Brain_Atlas_up <- dplyr::filter(listOfGO$Allen_Brain_Atlas_up,Adjusted.P.value <0.05)
  return(listOfGO)
}


goEnrichment <- lapply(res,dummyFun)
goEnrichmentFiltered <- lapply(goEnrichment,dummyFun2)

dataframeify <- function(df,colVal,nameOfCol){
  if(nrow(df)>0){
    df <- cbind(df,colVal,stringsAsFactors=F)
    colnames(df)[ncol(df)] <- nameOfCol
    df[,ncol(df)] 
  }else{
    df <- cbind(df,c())
    colnames(df)[ncol(df)] <- nameOfCol
  }
  return(df)
}
metaDFify <- function(a1){
  a2<-mapply(dataframeify,a1,names(a1),MoreArgs = list(nameOfCol='Library'),SIMPLIFY=F)
  a2 <- do.call(rbind,a2)
  return(a2)
}
a3 <- lapply(goEnrichmentFiltered,metaDFify)
a4 <- mapply(dataframeify,a3,names(a3),MoreArgs = list(nameOfCol='Branch'),SIMPLIFY=F)
a4 <- do.call(rbind,a4)


dv <- gsub('p_[0-9]\\.pos','UP',a4$Branch)
dv <- gsub('p_[0-9]\\.neg','DOWN',dv)
dlpfc_enrich <- dplyr::mutate(a4,Direction = dv)
mapFunction <- function(x){
  if(length(grep('1',x))>0){
    return(1)
  }else if (length(grep('2',x))>0){
    return(2)
  }else if (length(grep('3',x))>0){
    return(3)
  }else if (length(grep('4',x))>0){
    return(4)
  }else if (length(grep('5',x))>0){
    return(5)
  } else if (length(grep('6',x))>0){
    return(6)
  }
}
dlpfc_enrich <- dplyr::mutate(dlpfc_enrich,BranchNumber = sapply(dlpfc_enrich$Branch,mapFunction))
dlpfc_enrich <- dplyr::select(dlpfc_enrich,Term,Library,Direction,BranchNumber,P.value,Adjusted.P.value,Overlap,Combined.Score,Genes)
dlpfc_enrich$Overlap <- paste0('\'',dlpfc_enrich$Overlap)

#write.table(dlpfc_enrich,file='~/prot-lineage/results/F_GOenrichment.tsv',sep='\t',row.names=F,quote=F)
write.table(dlpfc_enrich,file='~/prot-lineage/results/M_GOenrichment.tsv',sep='\t',row.names=F,quote=F)






#for males there are 7 states
#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))
l2$p_6 <- rep(0,length(gene_short_name))
l2$p_7 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))
l2$d_7 <- rep(0,length(gene_short_name))


for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(MonRun$State2)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$p_3[i] <- tk$s[2,4]
  l2$p_4[i] <- tk$s[3,4]
  l2$p_5[i] <- tk$s[4,4]
  l2$p_6[i] <- tk$s[5,4]
  l2$p_7[i] <- tk$s[6,4]
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
  l2$d_6[i] <- tk$s[5,1]
  l2$d_7[i] <- tk$s[6,1]
}

#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))

dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)

dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])

df3 <- dplyr::left_join(dfa1,dfb1)

#gene names are currently labeled by gene and uniprot id. need to split:
df4 <- df3
df4$gene_short_name <- gsub("\\|.*", "", df4$gene_names)
#reorder columns
df4 <- df4[,c(1,5,2,3,4)]
#write.csv(df4,file='~/prot-lineage/results/female_anova_stats.csv',quote=F,row.names=F)
write.csv(df4,file='~/prot-lineage/results/male_anova_stats.csv',quote=F,row.names=F)





#run branch-specific enrichment analysis 
library(enrichR)

#use the initial data frame with all p-values & effect sizes from tukey analysis (df2)
dlpfc <- df2
dlpfc$gene_names <- gsub("\\|.*", "", dlpfc$gene_names)
getGeneSets = function(pval,dval,setName,setDescription,summaryMat){
  resPos <- dplyr::filter(summaryMat,.data[[pval]]<0.05, .data[[dval]]>0)
  resNeg <- dplyr::filter(summaryMat,.data[[pval]]<0.05, .data[[dval]]<0)
  res <- list()
  res$pos <- c(paste0(setName,'pos'),paste0(setDescription,'pos'),resPos$gene_names)
  res$neg <- c(paste0(setName,'neg'),paste0(setDescription,'neg'),resNeg$gene_names)
  return(res)
}
runList <- list()
runList$pval <- paste0('p_',2:7)
runList$dval <- paste0('d_',2:7)
runList$setName <- paste0('branch',2:7)
runList$setDescription <- paste0('branchSpecificDEG',2:7)
res <- mapply(getGeneSets,
              runList$pval,
              runList$dval,
              runList$setName,
              runList$setDescription,
              MoreArgs = list(summaryMat=dlpfc),
              SIMPLIFY= FALSE)
res <- unlist(res,recursive = F)

dummyFun <- function(geneSet){
  res <- enrichR::enrichr(geneSet[-c(1:2)],databases = c('GO_Biological_Process_2018',
                                                         'GO_Molecular_Function_2018',
                                                         'GO_Cellular_Component_2018'))
  return(res)
}
dummyFun2 <- function(listOfGO){
  listOfGO$GO_Biological_Process_2018 <- dplyr::filter(listOfGO$GO_Biological_Process_2018,Adjusted.P.value <0.05)
  listOfGO$GO_Molecular_Function_2018 <- dplyr::filter(listOfGO$GO_Molecular_Function_2018,Adjusted.P.value <0.05)
  listOfGO$GO_Cellular_Component_2018 <- dplyr::filter(listOfGO$GO_Cellular_Component_2018,Adjusted.P.value <0.05)
  return(listOfGO)
}




goEnrichment <- lapply(res,dummyFun)
goEnrichmentFiltered <- lapply(goEnrichment,dummyFun2)

dataframeify <- function(df,colVal,nameOfCol){
  if(nrow(df)>0){
    df <- cbind(df,colVal,stringsAsFactors=F)
    colnames(df)[ncol(df)] <- nameOfCol
    df[,ncol(df)] 
  }else{
    df <- cbind(df,c())
    colnames(df)[ncol(df)] <- nameOfCol
  }
  return(df)
}
metaDFify <- function(a1){
  a2<-mapply(dataframeify,a1,names(a1),MoreArgs = list(nameOfCol='Library'),SIMPLIFY=F)
  a2 <- do.call(rbind,a2)
  return(a2)
}
a3 <- lapply(goEnrichmentFiltered,metaDFify)
a4 <- mapply(dataframeify,a3,names(a3),MoreArgs = list(nameOfCol='Branch'),SIMPLIFY=F)
a4 <- do.call(rbind,a4)


dv <- gsub('p_[0-9]\\.pos','UP',a4$Branch)
dv <- gsub('p_[0-9]\\.neg','DOWN',dv)
dlpfc_enrich <- dplyr::mutate(a4,Direction = dv)
mapFunction <- function(x){
  if(length(grep('1',x))>0){
    return(1)
  }else if (length(grep('2',x))>0){
    return(2)
  }else if (length(grep('3',x))>0){
    return(3)
  }else if (length(grep('4',x))>0){
    return(4)
  }else if (length(grep('5',x))>0){
    return(5)
  } else if (length(grep('6',x))>0){
    return(6)
  } else if (length(grep('7',x))>0){
    return(7)
  }
}
dlpfc_enrich <- dplyr::mutate(dlpfc_enrich,BranchNumber = sapply(dlpfc_enrich$Branch,mapFunction))
dlpfc_enrich <- dplyr::select(dlpfc_enrich,Term,Library,Direction,BranchNumber,P.value,Adjusted.P.value,Overlap,Combined.Score,Genes)
dlpfc_enrich$Overlap <- paste0('\'',dlpfc_enrich$Overlap)

#write.table(dlpfc_enrich,file='~/prot-lineage/results/F_GOenrichment.tsv',sep='\t',row.names=F,quote=F)
write.table(dlpfc_enrich,file='~/prot-lineage/results/M_GOenrichment.tsv',sep='\t',row.names=F,quote=F)

plot_cell_trajectory(MonRun, color_by="cluster")
BEAM_res <- BEAM(MonRun, branch_point=1, cores=1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("fd","pval","qval")]
head(BEAM_res, n=10)
plot_genes_branched_heatmap(MonRun[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
MonRun_genes <- row.names(subset(fData(MonRun),
                               gene_short_name %in% c("CD63|F8VWK8", "LAMTOR4|C9JXA7", "CDC42SE2|Q9NRR3")))
plot_genes_branched_pseudotime(MonRun[MonRun_genes,],
                               branch_point = 1,
                               color_by = "diagnosis",
                               ncol = 1)
