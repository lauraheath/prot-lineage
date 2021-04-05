synapser::synLogin()
#Run script "get_proteins_script_JG.r" to get the log2 normalized matrix

dat <- Log2_Normalized
# dat <- sink_preWinzor
#dat <- Normalized
#dat <- 2^(Log2_Normalized)
dat[is.na(dat)] <- 0
dat <- as.matrix(dat)
dat <- Matrix(dat, sparse=TRUE)
gene_short_name <- as.data.frame(rownames(dat))
names(gene_short_name)[names(gene_short_name) == "rownames(dat)"] <- "gene_short_name"
rownames(gene_short_name) <- rownames(dat)

genesums <- data.frame(rowSums(dat), na.rm=T)
cellsums <- data.frame(colSums(dat), na.rm=T)
hist(genesums$rowSums.dat.)
hist(cellsums$colSums.dat.)

#impute missing pmi (2 missing: fill in with median pmi)
paste('Imputing PMI to:',median(Meta$pmi[!is.na(Meta$pmi)]))
#add this back into the metadata file
Meta$pmi[is.na(Meta$pmi)] <- 6.5


delog <- 2^(Log2_Normalized)
delog[is.na(delog)] <- 0
genesums <- data.frame(rowSums(delog), na.rm=T)
cellsums <- data.frame(colSums(delog), na.rm=T)
hist(genesums$rowSums.delog.)
hist(cellsums$colSums.delog.)

hist(delog$b01.129N)

#get list of genes that are differentially expressed between AD case & control (Jake's analysis)
#read in meta-analysis results (across 4 brain regions)
p <- synapser::synGet('syn22686552')
ADgenes_meta <- read.csv(p$path, sep="\t")
ADgenes_meta <- subset(ADgenes_meta, ADgenes_meta$fdr.random<0.05)
dim(ADgenes_meta)





ColNorm <- function(dat){

  M = max(colSums(dat))
  l <- length(colnames(dat))

  for( i in 1:l){

    dat[,i] = dat[,i]*(M/sum(dat[,i]))

  }

  return(dat)
}

datNorm <- ColNorm(dat)
cellsums2 <- data.frame(colSums(datNorm, na.rm=T))

dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- read.delim(dlpfcCPMObj$path,stringsAsFactors = F)
#to save the components of the monocle object for later use:
saveRDS(counts, file="LH_TCX_processed_countmatrix.rds")
saveRDS(gene_metadata, file="LH_TCX_processed_gene_metadata.rds")
saveRDS(temp2, file="LH_TCX_processed_cell_metadata.rds")
saveRDS(cds_tcx, file="LH_CDS_TCX_Monocle3Object.rds")

#for monocle: expression matrix=counts, cell_metadata=temp2, gene_metadata=gene_metadata

# must first unload synapser because causes multiple definitions of S4Vectors
detach("package:synapser", unload=TRUE)
unloadNamespace("PythonEmbedInR") 
#run monocle and get Monocle Object: cds_tcx


cds <- new_cell_data_set(datNorm, cell_metadata=Meta, gene_metadata=gene_short_name)

#limit to differentially expressed genes (via Jake's meta analysis across 4 regions)
genes2<-c()
for (gene in unique(c(as.vector(ADgenes_meta$peptide_id)))){
  if (gene %in% rownames(cds)){
    genes2 <- c(genes2,which(rownames(cds)==gene))
  }
}
length(genes2)
cds <- cds[genes2,]
dim(cds)

cds <- preprocess_cds(cds, method="PCA", norm_method="none", num_dim=30)
plot_pc_variance_explained(cds)
cds <- align_cds(cds, 
                 preprocess_method="PCA",
                 alignment_group="batch",
                 residual_model_formula_str="~educ+pmi")
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, cluster_method = "louvain")

cds <- learn_graph(cds)
#cds >- order_cells(cds)

cds$braaksc <- as.factor(cds$braaksc)
plot_cells(cds, color_cells_by="diagnosis", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
head(colData(cds))

plot_cells(cds, color_cells_by="batch", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds, color_cells_by="msex", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))


### create an indicator to root the tree nearest Control individuals with cogdx == 1


#separate by sex
cdsF <- cds[,cds$msex==0]
cdsM <- cds[,cds$msex==1]

cdsF <- clear_cds_slots(cdsF)
cdsF <- preprocess_cds(cdsF, method="PCA", norm_method="none", num_dim=30)
plot_pc_variance_explained(cdsF)
cdsF <- align_cds(cdsF, 
                 preprocess_method="PCA",
                 alignment_group="batch",
                 residual_model_formula_str="~educ+pmi")
cdsF <- reduce_dimension(cdsF, reduction_method = "UMAP")
cdsF <- cluster_cells(cdsF, cluster_method = "louvain")

cdsF <- learn_graph(cdsF)



#label clusters, root in the cluster with the most control patients & least AD patients
cdsF$cluster <- clusters(cdsF)
table(cdsF$cluster, cdsF$diagnosis)
prop.table(table(cdsF$cluster, cdsF$diagnosis))

#cluster 2 has highest proportion of control & low proportion of AD
plot_cells(cdsF, color_cells_by="diagnosis", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cdsF, color_cells_by="cluster", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#create an indicator for earliest/least pathology people (diagnosis==Control, ceradsc==4)
cdsF$early <- ifelse(cdsF$diagnosis=="Control" & cdsF$ceradsc==4, "early", "not_early")
get_earliest_principal_node <- function(cdsF, diagnosis="early"){
  cell_ids <- which(colData(cdsF)[, "early"]==diagnosis)
  
  closest_vertex <-
    cdsF@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cdsF), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cdsF)[["UMAP"]])$name[as.numeric(names
    (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
  
}
cdsF <- order_cells(cdsF, root_pr_nodes=get_earliest_principal_node(cdsF))
cdsF$pseudotime = pseudotime(cdsF)

cdsF$braaksc <- as.factor(cdsF$braaksc)
cdsF$ceradsc <- as.factor(cdsF$ceradsc)
cdsF$cogdx <- as.factor(cdsF$cogdx)
plot_cells(cdsF, color_cells_by="pseudotime", cell_size = 2)
plot_cells(cdsF, color_cells_by="diagnosis", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
head(colData(cdsF))
plot_cells(cdsF, color_cells_by="braaksc", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cdsF, color_cells_by="cogdx", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cdsF, color_cells_by="ceradsc", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
plot_cells(cdsF, color_cells_by="batch", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#Make new dataframe with SampleID,Pseudotime, sex, diagnosis, braak, ceradsc, cogdx, apoe, pmi, age_death
Fvariables <- list()
Fvariables$SampleID <- cdsF$SampleID
Fvariables$pseudotime <- cdsF$pseudotime
Fvariables$msex <- cdsF$msex
Fvariables$diagnosis <- cdsF$diagnosis
Fvariables$braaksc <- cdsF$braaksc
Fvariables$ceradsc <- cdsF$ceradsc
Fvariables$cogdx <- cdsF$cogdx
Fvariables$apoe_genotype <- cdsF$apoe_genotype
Fvariables$APOE <- cdsF$APOE
Fvariables$pmi <- cdsF$pmi
Fvariables$age_death <- cdsF$age_death
Fvariables <- as.data.frame(Fvariables)
Fvariables$pseudotime_sc <- scale(Fvariables$pseudotime, center=F)

casecontrolF <- subset(Fvariables, Fvariables$diagnosis=='AD'|Fvariables$diagnosis=='Control')
casecontrolF$diag2 <- ifelse(casecontrolF$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ pseudotime_sc,casecontrolF,family='binomial'))
ggplot(casecontrolF,aes(x=diagnosis,
                      y=pseudotime_sc,
                      color=diagnosis)) + geom_boxplot()


hist(Fvariables$pseudotime)
table(Fvariables$diagnosis, Fvariables$braaksc)
table(Fvariables$diagnosis, Fvariables$cogdx)
table(Fvariables$diagnosis, Fvariables$ceradsc)
table(Fvariables$cogdx, Fvariables$ceradsc)
table(Fvariables$cogdx, Fvariables$braaksc)


braakfit <- MASS::polr(braaksc ~ pseudotime,Fvariables)
ceradfit <- MASS::polr(ceradsc ~ pseudotime,Fvariables)
cogdxfit <- MASS::polr(cogdx ~ pseudotime,Fvariables)

cat('braak p-value: ',pt(abs(summary(braakfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
cat('cerad p-value: ',pt(abs(summary(ceradfit)$coef[1,3]),ceradfit$df.residual,lower.tail=F)*2,'\n')
cat('cogdx p-value: ',pt(abs(summary(cogdxfit)$coef[1,3]),cogdxfit$df.residual,lower.tail=F)*2,'\n')


g <- ggplot2::ggplot(Fvariables, aes(x=braaksc, y=scale(pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g

g <- ggplot2::ggplot(Fvariables, aes(x=ceradsc, y=scale(pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g

g <- ggplot2::ggplot(Fvariables, aes(x=cogdx, y=scale(pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="COGDX",y="Pseudotime",x="COGDX")
g





######## male samples ############

cdsM <- clear_cds_slots(cdsM)
cdsM <- preprocess_cds(cdsM, method="PCA", norm_method="none", num_dim=30)
plot_pc_variance_explained(cdsM)
cdsM <- align_cds(cdsM, 
                 preprocess_method="PCA",
                 alignment_group="batch",
                 residual_model_formula_str="~educ+pmi")
cdsF <- reduce_dimension(cdsF, reduction_method = "UMAP")
cdsM <- reduce_dimension(cdsM, reduction_method = "UMAP")
cdsM <- cluster_cells(cdsM, cluster_method = "louvain")

cdsM <- learn_graph(cdsM)
#create an indicator for earliest/least pathology people (diagnosis==Control, ceradsc==4)
cdsM$braaksc <- as.numeric(cdsM$braaksc)
cdsM$early <- ifelse(cdsM$diagnosis=="Control" & cdsM$ceradsc==4 & cdsM$braaksc<3, "early", "not_early")
get_earliest_principal_node <- function(cdsM, diagnosis="early"){
  cell_ids <- which(colData(cdsM)[, "early"]==diagnosis)
  
  closest_vertex <-
    cdsM@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cdsM), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cdsM)[["UMAP"]])$name[as.numeric(names
                                                               (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
  
}

cdsM <- order_cells(cdsM, root_pr_nodes=get_earliest_principal_node(cdsM))
cdsM$pseudotime = pseudotime(cdsM)

cdsM$braaksc <- as.factor(cdsM$braaksc)
cdsM$ceradsc <- as.factor(cdsM$ceradsc)
cdsM$cogdx <- as.factor(cdsM$cogdx)
plot_cells(cdsM, color_cells_by="pseudotime", cell_size = 2)
plot_cells(cdsM, color_cells_by="diagnosis", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
head(colData(cdsM))
plot_cells(cdsM, color_cells_by="braaksc", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cdsM, color_cells_by="cogdx", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#Make new dataframe with SampleID,Pseudotime, sex, diagnosis, braak, ceradsc, cogdx, apoe, pmi, age_death
Mvariables <- list()
Mvariables$SampleID <- cdsM$SampleID
Mvariables$pseudotime <- cdsM$pseudotime
Mvariables$msex <- cdsM$msex
Mvariables$diagnosis <- cdsM$diagnosis
Mvariables$braaksc <- cdsM$braaksc
Mvariables$ceradsc <- cdsM$ceradsc
Mvariables$cogdx <- cdsM$cogdx
Mvariables$apoe_genotype <- cdsM$apoe_genotype
Mvariables$APOE <- cdsM$APOE
Mvariables$pmi <- cdsM$pmi
Mvariables$age_death <- cdsM$age_death
Mvariables <- as.data.frame(Mvariables)
Mvariables$pseudotime_sc <- scale(Mvariables$pseudotime, center=F)

casecontrolM <- subset(Mvariables, Mvariables$diagnosis=='AD'|Mvariables$diagnosis=='Control')
casecontrolM$diag2 <- ifelse(casecontrolM$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ pseudotime_sc,casecontrolM,family='binomial'))
ggplot(casecontrolM,aes(x=diagnosis,
                        y=pseudotime_sc,
                        color=diagnosis)) + geom_boxplot()


hist(Mvariables$pseudotime)
table(Mvariables$diagnosis, Mvariables$braaksc)
table(Mvariables$diagnosis, Mvariables$cogdx)
table(Mvariables$diagnosis, Mvariables$ceradsc)
table(Mvariables$cogdx, Mvariables$ceradsc)
table(Mvariables$cogdx, Mvariables$braaksc)


braakfit <- MASS::polr(braaksc ~ pseudotime,Mvariables)
ceradfit <- MASS::polr(ceradsc ~ pseudotime,Mvariables)
cogdxfit <- MASS::polr(cogdx ~ pseudotime,Mvariables)

cat('braak p-value: ',pt(abs(summary(braakfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
cat('cerad p-value: ',pt(abs(summary(ceradfit)$coef[1,3]),ceradfit$df.residual,lower.tail=F)*2,'\n')
cat('cogdx p-value: ',pt(abs(summary(cogdxfit)$coef[1,3]),cogdxfit$df.residual,lower.tail=F)*2,'\n')


g <- ggplot2::ggplot(Mvariables, aes(x=braaksc, y=scale(pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g

g <- ggplot2::ggplot(Mvariables, aes(x=ceradsc, y=scale(pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g

g <- ggplot2::ggplot(Mvariables, aes(x=cogdx, y=scale(pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="COGDX",y="Pseudotime",x="COGDX")
g




######## ALL samples ############

cds <- clear_cds_slots(cds)
cds <- preprocess_cds(cds, method="PCA", norm_method="none", num_dim=30)
plot_pc_variance_explained(cds)
cds <- align_cds(cds, 
                 preprocess_method="PCA",
                 alignment_group="batch",
                 residual_model_formula_str="~educ+pmi")
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, cluster_method = "louvain")
dim(cds)
cds <- learn_graph(cds)
#create an indicator for earliest/least pathology people (diagnosis==Control, ceradsc==4)

cds$early <- ifelse(cds$diagnosis=="Control" & cds$ceradsc==4, "early", "not_early")
get_earliest_principal_node <- function(cds, diagnosis="Control"){
  cell_ids <- which(colData(cds)[, "diagnosis"]==diagnosis)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
  
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
cds$pseudotime = pseudotime(cds)

cds$braaksc <- as.factor(cds$braaksc)
cds$ceradsc <- as.factor(cds$ceradsc)
cds$cogdx <- as.factor(cds$cogdx)
plot_cells(cds, color_cells_by="pseudotime", cell_size = 2)
plot_cells(cds, color_cells_by="diagnosis", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
head(colData(cds))
plot_cells(cds, color_cells_by="braaksc", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds, color_cells_by="cogdx", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

plot_cells(cds, color_cells_by="ceradsc", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))

#Make new dataframe with SampleID,Pseudotime, sex, diagnosis, braak, ceradsc, cogdx, apoe, pmi, age_death
variables <- list()
variables$SampleID <- cds$SampleID
variables$pseudotime <- cds$pseudotime
variables$msex <- cds$msex
variables$diagnosis <- cds$diagnosis
variables$braaksc <- cds$braaksc
variables$ceradsc <- cds$ceradsc
variables$cogdx <- cds$cogdx
variables$apoe_genotype <- cds$apoe_genotype
variables$APOE <- cds$APOE
variables$pmi <- cds$pmi
variables$age_death <- cds$age_death
variables <- as.data.frame(variables)
variables$pseudotime_sc <- scale(variables$pseudotime, center=F)

casecontrol <- subset(variables, variables$diagnosis=='AD'|variables$diagnosis=='Control')
casecontrol$diag2 <- ifelse(casecontrol$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ pseudotime_sc,casecontrol,family='binomial'))
ggplot(casecontrol,aes(x=diagnosis,
                       y=pseudotime_sc,
                       color=diagnosis)) + geom_boxplot()


hist(variables$pseudotime)
table(variables$diagnosis, variables$braaksc)
table(variables$diagnosis, variables$cogdx)
table(variables$diagnosis, variables$ceradsc)
table(variables$cogdx, variables$ceradsc)
table(variables$cogdx, variables$braaksc)


braakfit <- MASS::polr(braaksc ~ pseudotime,variables)
ceradfit <- MASS::polr(ceradsc ~ pseudotime,variables)
cogdxfit <- MASS::polr(cogdx ~ pseudotime,variables)

cat('braak p-value: ',pt(abs(summary(braakfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
cat('cerad p-value: ',pt(abs(summary(ceradfit)$coef[1,3]),ceradfit$df.residual,lower.tail=F)*2,'\n')
cat('cogdx p-value: ',pt(abs(summary(cogdxfit)$coef[1,3]),cogdxfit$df.residual,lower.tail=F)*2,'\n')


g <- ggplot2::ggplot(variables, aes(x=braaksc, y=scale(pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g

g <- ggplot2::ggplot(variables, aes(x=ceradsc, y=scale(pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g

g <- ggplot2::ggplot(variables, aes(x=cogdx, y=scale(pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="COGDX",y="Pseudotime",x="COGDX")
g
