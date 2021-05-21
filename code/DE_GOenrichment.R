######### branch-specific differential expression analysis

#need Monrun object from protein_lineage_monocle2 and gene_short_name vector
#script for males following script for females (7 states in males, 6 in females)

MonRun <- readRDS(file="~/prot-lineage/results/Female_monocleObject.rds")
temp <- readRDS(file="~/prot-lineage/results/Female_matrix.rds")
#there are 6 states in the female tree
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

write.csv(df4,file='~/prot-lineage/results/female_DEanova_stats.csv',quote=F,row.names=F)


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
  }
}
dlpfc_enrich <- dplyr::mutate(dlpfc_enrich,BranchNumber = sapply(dlpfc_enrich$Branch,mapFunction))
dlpfc_enrich <- dplyr::select(dlpfc_enrich,Term,Library,Direction,BranchNumber,P.value,Adjusted.P.value,Overlap,Combined.Score,Genes)
dlpfc_enrich$Overlap <- paste0('\'',dlpfc_enrich$Overlap)

write.table(dlpfc_enrich,file='~/prot-lineage/results/F_GOenrichment.tsv',sep='\t',row.names=F,quote=F)





#run Male monocle object with script below:
MonRun <- readRDS(file="~/prot-lineage/results/Male_monocleObject.rds")
temp <- readRDS(file="~/prot-lineage/results/Male_matrix.rds")

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

write.csv(df4,file='~/prot-lineage/results/male_DEanova_stats.csv',quote=F,row.names=F)

#run branch-specific enrichment analysis 

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

write.table(dlpfc_enrich,file='~/prot-lineage/results/M_GOenrichment.tsv',sep='\t',row.names=F,quote=F)

