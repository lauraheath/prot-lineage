synapser::synLogin()
#Run script "get_proteins_script_JG.r" to get the log2 normalized matrix

dat <- Log2_Normalized
dat <- sink_preWinzor
dat <- Normalized
dat[is.na(dat)] <- 0
dat <- as.matrix(dat)
dat <- Matrix(dat, sparse=TRUE)
gene_short_name <- as.data.frame(rownames(dat))
names(gene_short_name)[names(gene_short_name) == "rownames(dat)"] <- "gene_short_name"
rownames(gene_short_name) <- rownames(dat)

genesums <- data.frame(rowSums(dat))
cellsums <- data.frame(colSums(dat))

dat <- data.matrix(dat)
hist(dat)

ColNorm <- function(dat){

  M = max(colSums(dat))
  l <- length(colnames(dat))

  for( i in 1:l){

    dat[,i] = dat[,i]*(M/sum(dat[,i]))

  }

  return(dat)
}

datNorm <- ColNorm(dat)
cellsums2 <- data.frame(colSums(datNorm))
hist(datNorm)

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

cds <- preprocess_cds(cds, method="PCA", norm_method="log", num_dim=50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, cluster_method = "louvain")

cds <- learn_graph(cds)
#cds >- order_cells(cds)

cds$braaksc <- as.factor(cds$braaksc)
plot_cells(cds, xcolor_cells_by="pseudotime")
plot_cells(cds, color_cells_by="diagnosis", cell_size = 2, label_cell_groups=0,show_trajectory_graph=TRUE)+theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 7),legend.key.size = unit(.5, "cm"),
  legend.key.width = unit(0.5,"cm"))
head(colData(cds))
 #separate by sex
cdsF <- cds[,cds$msex==0]
cdsM <- cds[,cds$msex==1]

cdsF <- clear_cds_slots(cdsF)
cdsF <- preprocess_cds(cdsF, method="PCA", norm_method="none", num_dim=50)
plot_pc_variance_explained(cdsF)
cdsF <- reduce_dimension(cdsF, reduction_method = "UMAP")
cdsF <- cluster_cells(cdsF, cluster_method = "louvain")

cdsF <- learn_graph(cdsF)
get_earliest_principal_node <- function(cdsF, Diagnosis="Control"){
  cell_ids <- which(colData(cdsF)[, "diagnosis"]==Diagnosis)
  
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

casecontrolF <- subset(Fvariables, Fvariables$diagnosis=='AD'|Fvariables$diagnosis=='Control')
casecontrolF$diag2 <- ifelse(casecontrolF$diagnosis=='AD', 1, 0)


casecontrolF$pseudotime_sc <- scale(casecontrolF$pseudotime,center=F)
summary(glm(diag2 ~ pseudotime_sc,casecontrolF,family='binomial'))
ggplot(casecontrolF,aes(x=diagnosis,
                      y=pseudotime,
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
cdsM <- preprocess_cds(cdsM, method="PCA", norm_method="none", num_dim=50)
plot_pc_variance_explained(cdsM)
cdsM <- reduce_dimension(cdsM, reduction_method = "UMAP")
cdsM <- cluster_cells(cdsM, cluster_method = "louvain")

cdsM <- learn_graph(cdsM)
get_earliest_principal_node <- function(cdsM, Diagnosis="Control"){
  cell_ids <- which(colData(cdsM)[, "diagnosis"]==Diagnosis)
  
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
x <- list()
x$SampleID <- cdsM$SampleID
x$pseudotime <- cdsM$pseudotime
x$msex <- cdsM$msex
x$diagnosis <- cdsM$diagnosis
x$braaksc <- cdsM$braaksc
x$ceradsc <- cdsM$ceradsc
x$cogdx <- cdsM$cogdx
x$apoe_genotype <- cdsM$apoe_genotype
x$APOE <- cdsM$APOE
x$pmi <- cdsM$pmi
x$age_death <- cdsM$age_death
x <- as.data.frame(x)

casecontrolM <- subset(x, x$diagnosis=='AD'|x$diagnosis=='control')

casecontrolM$pseudotime_sc <- scale(casecontrolM$pseudotime,center=F)
summary(glm(diagnosis ~ pseudotime_sc,casecontrolM,family='binomial'))
summary(glm(diagnosis ~ pseudotime,casecontrolM,family='binomial'))


malefemale <- rbind(casecontrolF, casecontrolM)
malefemale$Sex <- ifelse(malefemale$msex==0, "Female", "Male")
ggplot(malefemale,aes(x=Sex,
                      y=pseudotime,
                      color=diagnosis)) + geom_boxplot()


