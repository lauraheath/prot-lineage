

#get bulk rnaseq pseudotimes and see if proteomic-derived pseudotimes are correlated with rna-seq pseudotimes, for those patients
#included in both studies (n=107 in females)
#need patient IDs across studies (run "get_proteins_metadata.r" to generate this file in your local directory)
patients <- read.csv(file="~/prot-lineage/data/patientIDS_rnaseq_proteomics.csv")

#first run protein_lineage_monocle2.R to get covariates file, or read in variables file from data folder if saved
Fvariables <- read.csv(file="~/prot-lineage/results/prot_pstime_covars_F.csv")
MonRun <- readRDS(file="~/prot-lineage/results/Female_monocleObject.rds")
#males
#Fvariables <- read.csv(file="~/prot-lineage/results/prot_pstime_covars_M.csv")
#MonRun <- readRDS(file="~/prot-lineage/results/Male_monocleObject.rds")


#also run rerun_rnaseq_lineage.r or use saved gene list from that script:
rnaseq_genes <- read.csv(file="~/prot-lineage/data/rnaseq_genes.csv")
prot_genes <- read.csv(file="~/prot-lineage/data/prot_genes.csv")

#venn diagram of rnaseq genes and prot genes
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")


rna <- rnaseq_genes$gene_short_name
rna <- rnaseq_genes[!(rnaseq_genes$gene_short_name==""),]
prot_genes$gene_name <- gsub("\\|.*", "", prot_genes$gene_short_name)
prot <- prot_genes$gene_name


venn.diagram(
  x=list(rna, prot),
  category.names=c("RNAseq Genes","Proteins"),
  filename='Gene_protein_overlap_venn.png',
  output=TRUE,
  imagetype="png",
  height=1500,
  width=2000,
  resolution=300,
  compression="lzw",
  lwd = 2,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.default.post="outer",
  cat.pos=c(-20,20),
  cat.dist=c(0.05,0.05)
)


#venn diagram of overlapping patients in rnaseq and protein studies

patients$individualID<-as.character(patients$individualID)
rna_ps <- subset(patients, patients$rnaseq==1)
prot_ps <- subset(patients, patients$proteomics==1)
rna_patients <- rna_ps$individualID
prot_patients <- prot_ps$individualID

venn.diagram(
  x=list(rna_patients, prot_patients),
  category.names=c("RNAseq IDs","Proteomics IDs"),
  filename='Patient_overlap_venn.png',
  output=TRUE,
  imagetype="png",
  height=1500,
  width=2000,
  resolution=300,
  compression="lzw",
  lwd = 2,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.default.post="outer",
  cat.pos=c(-24,18),
  cat.dist=c(0.05,0.05)
)





Fvariables$braaksc <- factor(Fvariables$braaksc,levels = c(0:6))
Fvariables$ceradsc <- factor(Fvariables$ceradsc,levels = c(1:4))
Fvariables$cogdx <- factor(Fvariables$cogdx,levels = c(1:6))

#upload pseudotimes from bulk rna-seq lineage paper
#females
p <- synapser::synGet('syn22822695')
rnaseq_pstimes <- read.csv(p$path)
#males
#p <- synapser::synGet('syn23446654')
#rnaseq_pstimes <- read.csv(p$path)

rnaseq_pstimes$rnaseq_pseudotime_sc <- scale(rnaseq_pstimes$Pseudotime, center=F)
names(rnaseq_pstimes)[names(rnaseq_pstimes) == "SampleID"] <- "rnaseqID"
names(rnaseq_pstimes)[names(rnaseq_pstimes) == "State"] <- "rnaseq_State"
names(rnaseq_pstimes)[names(rnaseq_pstimes) == "Pseudotime"] <- "rnaseq_Pseudotime"
rnaseq_pstimes$roROW_ID<-NULL
rnaseq_pstimes$ROW_VERSION<-NULL
rnaseq_pstimes$col<-NULL

corrs <- merge(rnaseq_pstimes, Fvariables, by="rnaseqID")
cor(corrs$rnaseq_pseudotime_sc, corrs$pseudotime_sc, method="pearson")
cor(corrs$rnaseq_pseudotime_sc, corrs$pseudotime_sc, method="spearman")
plot(corrs$rnaseq_pseudotime_sc, corrs$pseudotime_sc)
abline(lm(corrs$rnaseq_pseudotime_sc~corrs$pseudotime_sc))
lines(lowess(corrs$rnaseq_pseudotime_sc,corrs$pseudotime_sc))
corrs$braaksc<-as.factor(corrs$braaksc)
corrs$ceradsc<-as.factor(corrs$ceradsc)
corrs$cogdx<-as.factor(corrs$cogdx)

#add rna-seq pseudotimes to the monocle object
rnaseq_pstime <- subset(corrs, select=c(SampleID, rnaseq_pseudotime_sc))
Fvariables2 <- merge(Fvariables, rnaseq_pstime, all=TRUE)
MonRun$rnaseq_pseudotime <- Fvariables2$rnaseq_pseudotime_sc
MonRun$prot_pseudotime <- Fvariables2$pseudotime_sc


#tiff(file='~/prot-lineage/figures/FEMALE_tree_rnaseq_pstime.tiff',height=150,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_tree_rnaseq_pstime.tiff',height=150,width=100,units='mm',res=300)
g <- plot_cell_trajectory(MonRun,color_by = "rnaseq_pseudotime",show_branch_points=F,use_color_gradient = F,cell_size = 1) + scale_color_gradient (low="blue", high="red")
g
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="pseudotime_sc", y="rnaseq_pseudotime_sc", 
                  add="reg.line",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_diagnosis.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_diagnosis.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="pseudotime_sc", y="rnaseq_pseudotime_sc", color = "diagnosis",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_braak.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="pseudotime_sc", y="rnaseq_pseudotime_sc", color = "braaksc",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_cerad.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="pseudotime_sc", y="rnaseq_pseudotime_sc", color = "ceradsc",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_cogdx.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_cogdx.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="pseudotime_sc", y="rnaseq_pseudotime_sc", color = "cogdx",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()





#compare DE genes included as feature sets in bulk rnaseq and prot

prot_genes$gene_short_name <- gsub("\\|.*", "", prot_genes$gene_short_name)
inters <- intersect(rnaseq_genes$gene_short_name,prot_genes$gene_short_name)
diffs <- setdiff(rnaseq_genes$gene_short_name,prot_genes$gene_short_name)

#rerun lineage analysis on both rna-seq and proteomics data using only the shared genes as a feature set
inters <- as.data.frame(inters)
names(inters)[names(inters) == "inters"] <- "gene_short_name"
Log2_Normalized <- readRDS(file="~/prot-lineage/data/Log2_Normalized.rds")
Meta <- readRDS(file="~/prot-lineage/data/Meta.rds")
Log2_Normalized2 <- Log2_Normalized
Log2_Normalized2$proteins <- rownames(Log2_Normalized2)
Log2_Normalized2$proteins <- gsub("\\|.*", "", Log2_Normalized2$proteins)

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

Dat <- Log2_Normalized2
Dat[is.na(Dat)] <- 0

#select only the rows with gene short names that match the intersecting gene short names (some are repeated peptides, same gene)
genes2<-c()
for (gene in unique(c(as.vector(inters$gene_short_name)))){
  if (gene %in% Dat$proteins){
    genes2 <- c(genes2,which(Dat$proteins==gene))
  }
}
length(genes2)
Dat2 <- Dat[genes2,]
dim(Dat2)
Dat2$proteins<-NULL


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

MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g

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
F_inters_prots <- as.data.frame(x)
F_inters_prots$pseudotime_sc <- scale(F_inters_prots$Pseudotime, center=F)




#### now upload rna-seq data and rerun lineage analysis on intersecting proteins


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

dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- read.delim(dlpfcCPMObj$path,stringsAsFactors = F)

#synapse id of dat2 file (rosmap covariates): syn8466814
dlpfcCovObj <- synapser::synGet('syn11024258')
covars <- read.delim(dlpfcCovObj$path,stringsAsFactors = F)

#converting ENSG to gene symbols
convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           #host='useast.ensembl.org')
                           host='uswest.ensembl.org')
  
  genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               filters='ensembl_gene_id',
               values=ensemblIds,
               mart=ensembl)
  return(genes)
}
Make.Gene.Symb <- function(GeneENSG){
  
  #source('convertEnsemblToHgnc.R')
  GeneConv <- convertEnsemblToHgnc(GeneENSG)
  Symb <- as.character(c(1:length(GeneENSG)))
  
  for (i in 1:length(GeneENSG)){
    In <- which(GeneConv$ensembl_gene_id == GeneENSG[i])
    if (length(In)>0){
      Symb[i] <- GeneConv$external_gene_name[In]
    }
  }
  
  return(Symb)
  
}
Dat$gene_short_name <- Make.Gene.Symb(Dat$ensembl_gene_id)

#select only the rows with gene short names that match the intersecting gene short names (some are repeated peptides, same gene)
genes2<-c()
for (gene in unique(c(as.vector(inters$gene_short_name)))){
  if (gene %in% Dat$gene_short_name){
    genes2 <- c(genes2,which(Dat$gene_short_name==gene))
  }
}
length(genes2)
Dat2 <- Dat[genes2,]
dim(Dat2)


Names <- colnames(Dat2)

for (i in 1:length(Names)){
  
  Names[i] <- substring(Names[i],2)
  
}


colnames(Dat2) <- Names
cNames <- covars$SampleID
l <- length(Names)

#deleting columns not in the covariate list
temp <- rep(T,l)
for (i in 1:l){
  if (!(Names[i] %in% cNames)){
    temp[i] <- F
  }
}

In <- which(temp)
#print(temp)
Dat2 <- Dat2[,In]

#deleting extra rows in covariate list
Names <- Names[In]
l <- length(cNames)
temp <- rep(T,l)
for (i in 1:l){
  if (!(cNames[i] %in% Names)){
    temp[i] <- F
  }
}
In <- which(temp)
covars <- covars[In,]

ColNorm <- function(Dat2){
  
  M = max(colSums(Dat2))
  l <- length(colnames(Dat2))
  
  for( i in 1:l){
    
    Dat2[,i] = Dat2[,i]*(M/sum(Dat2[,i]))
    
  }
  
  return(Dat2)
}

DatNorm <- ColNorm(Dat2)

#removing bad batches
DatNorm <- DatNorm[,covars$Batch<7]
covars <- covars[covars$Batch<7,] 


#Keeping only female data 
#Sex <- 'FEMALE'
In_S <- which(covars$msex == 0)
DatNorm2 <- DatNorm[,In_S]
covars2 <- covars[In_S,]

temp <- DatNorm2
temp2 <- covars2

rosmapObj <- synapser::synGet('syn3191087')
rosmap <- data.table::fread(rosmapObj$path,data.table=F)

#add in braak score & cerad score
rosmapIdObj <- synapser::synGet('syn3382527')
rosmapId <- data.table::fread(rosmapIdObj$path,data.table=F)
rosmapId <- dplyr::select(rosmapId,projid,rnaseq_id)
rosmapRNAid<-dplyr::left_join(rosmapId,rosmap)
#remove duplicate rows
rosmapRNAid <- unique(rosmapRNAid)
rosmapRNAid2 <- subset(rosmapRNAid, select=c(rnaseq_id,braaksc,ceradsc))
names(rosmapRNAid2)[names(rosmapRNAid2) == "rnaseq_id"] <- "SampleID"

temp2<-dplyr::left_join(temp2,rosmapRNAid2, by="SampleID")

temp2$braaksc <- factor(temp2$braaksc,levels = c(0:6))
temp2$ceradsc <- factor(temp2$ceradsc,levels = c(1:4))
temp2$cogdx <- factor(temp2$cogdx, levels = c(1:6))


gene_short_name <- inters$gene_short_name

rownames(temp)<-NULL
rownames(temp2)<-NULL

MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g

x <- list()
x$rnaseqID <- MonRun$SampleID
x$State_rnaseq <- MonRun$State
x$Pseudotime_rnaseq <- MonRun$Pseudotime
F_inters_rnaseq <- as.data.frame(x)
F_inters_rnaseq$pseudotime_scRNA <- scale(F_inters_rnaseq$Pseudotime, center=F)


F_intersecting_pseudotimes <- merge(F_inters_prots, F_inters_rnaseq, by="rnaseqID")

cor(F_intersecting_pseudotimes$pseudotime_sc, F_intersecting_pseudotimes$pseudotime_scRNA, method="pearson")
cor(F_intersecting_pseudotimes$pseudotime_sc, F_intersecting_pseudotimes$pseudotime_scRNA, method="spearman")
plot(F_intersecting_pseudotimes$pseudotime_sc, F_intersecting_pseudotimes$pseudotime_scRNA)
abline(lm(F_intersecting_pseudotimes$pseudotime_sc~F_intersecting_pseudotimes$pseudotime_scRNA))
lines(lowess(corrs$rnaseq_pseudotime_sc,corrs$pseudotime_sc))

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="Pseudotime", y="rnaseq_Pseudotime", 
                  add="reg.line",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()
