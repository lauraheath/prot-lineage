devtools::install_github("kassambara/ggpubr")
library(ggpubr)

#load data functions
#' SynID <- 'syn23583548'
#' A function
#' 
#' @param SynID A vector
#' @param names A numeric
#' 
#' @return The dataframe pulled from synapse
#' @importFrom synapser synGet
#' @importFrom data.table fread
#' @export

TMT_Express_Load <- function( SynID, names ){
  #'@param  SynID a synapse ID of a proteomics csv matrix eg. syn21266454 or syn21266454
  #'@param  names number corresponding to the row that are the rownames (Use 0 if none apply )
  if( !is.character(SynID) ){
    return("Invalid Synapse ID: Input SynID is not a character string")
  }
  if( !grepl('syn', SynID) | !((nchar(SynID) == 11) | (nchar(SynID) == 10 ))  ){
    return("Invalid Synapse ID: Input SynID is not a valid synapse ID")
  }
  #if(  grepl('Benefactor not found for', as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) )) ){
  #  return("Syn ID does not exist")
  #}
  if(  'You do not have READ permission for the requested entity' %in% as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) ) ){
    return("User does not have READ access")
  }
  if( !grepl('.csv$', synapser::synGet(SynID)$path ) ){
    return("File to import must be a CSV File")
  }
  if( !is.numeric(names) ){
    return("names is not a number")
  }
  
  if( names == 0 ){
    import <- data.frame( data.table::fread( synapser::synGet(SynID)$path, header=T, sep=',' ))
  }else{
    import <- data.frame( data.table::fread( synapser::synGet(SynID)$path, header=T, sep=',' ), row.names = 1)
  }
  return( import )
}




## Import Data
Log2_Normalized <-TMT_Express_Load('syn21266454', 1)
Normalized <- TMT_Express_Load('syn21266453', 1)
Meta <- TMT_Express_Load('syn21323404', 1)
# - Public Facing BioSpecimin Data: syn21323366
# - Staged BioSpecimin Data: syn23583548
BioSpecimin <- TMT_Express_Load('syn21323366', 0)
if( "assay" %in% colnames(BioSpecimin) ){
}else{
  BioSpecimin <- TMT_Express_Load('syn23583548', 0)
}
BioSpecimin <- BioSpecimin[ BioSpecimin$assay == 'TMT quantitation', ]

Meta$specimenID <- row.names(Meta)
Meta <- dplyr::left_join(Meta, BioSpecimin, by = 'specimenID' )

Clinical <- TMT_Express_Load( 'syn3191087',0 )
Meta <- dplyr::left_join(Meta, Clinical, by = 'individualID' )
Meta <- Meta[ ,colSums(is.na(Meta))<nrow(Meta) ] 
#remove pool samples
Meta <- Meta[ Meta$isAssayControl==FALSE,]
row.names( Meta ) <- Meta$batchChannel
Meta <- Meta[ colnames(Log2_Normalized), ]
Meta <- Meta[ ,colnames(Meta)[ (colnames(Meta) %in%'controlType' )==F] ]
# Harmonize case-control status
Meta$braaksc <- as.numeric(Meta$braaksc)
Meta$ceradsc <- as.numeric(Meta$ceradsc)
Meta$cogdx <- as.numeric(Meta$cogdx)
# Harmonize case-control status
Meta$diagnosis <- "other"
Meta$diagnosis[Meta$cogdx == 1 & Meta$braaksc <= 3 & Meta$ceradsc >= 3] <- "control"
Meta$diagnosis[Meta$cogdx == 4 & Meta$braaksc >= 4 & Meta$ceradsc <= 2] <- "AD"
kableExtra::kable( table(Meta$diagnosis) )
## Add Ages over 90 for modeling
Mast <- TMT_Express_Load('syn23573928', 0)
Mast<- Mast[ !duplicated(Mast$individualID), ]
Meta <- dplyr::left_join(Meta[ , colnames(Meta)[ (colnames(Meta) %in% c('age_at_visit_max', 'age_first_ad_dx', 'age_death') )==F ] ], Mast[, c('individualID', 'age_at_visit_max', 'age_first_ad_dx', 'age_death')],by = 'individualID')
#Convert APOE
Meta$apoe_genotype <- as.numeric( Meta$apoe_genotype )
APOS <- as.numeric(names( table(Meta$apoe_genotype) ))
names(APOS) <- APOS
APOS[names( table(Meta$apoe_genotype) )] <- 0
APOS[grepl("4",names(APOS))] <- 1
APOS[grepl("44",names(APOS))] <- 2
Meta$APOE <- as.numeric( APOS[ as.character(Meta$apoe_genotype) ] )
## Winzorize Expression Data
sink_preWinzor<- Log2_Normalized
for( i in 1:dim(Log2_Normalized)[1] ){
  Log2_Normalized[i,] <- DescTools::Winsorize( as.numeric(Log2_Normalized[i,]), na.rm = TRUE ) 
}
row.names(Meta) <- Meta$batchChannel
#Code Sex and Diagnosis as Factors
Meta$msex <- as.factor(Meta$msex)
Meta$APOE <- as.factor(Meta$APOE)
####Meta for Diagnosis
Meta_D <- Meta[ Meta$diagnosis %in% c('AD','control'), ]
Meta_D$diagnosis <- as.factor(Meta_D$diagnosis)


#pare log2 matrix to AD case and controls
ADmatrix <- Log2_Normalized[Meta_D$batchChannel]

#run limma 
design <- model.matrix(~ 0 + diagnosis + msex + pmi, data = Meta_D)
design
cont.matrix <- makeContrasts(diagnosisAD-diagnosiscontrol, levels=colnames(design))
cont.matrix
fit <- lmFit(ADmatrix, design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)
topTable(fit)

#extract results for all genes and get 95% confidence limits
allresults <- topTable(fit, n=Inf, confint=0.95)

#relabel columns to match lfq data for agora
#Fix protein names with most updated
p1 <- synapser::synGet('syn24216770')
correct_geneIDs <- read.csv(p1$path)
allresults$OldPeptideID <- rownames(allresults)
allresults <- dplyr::left_join(allresults, correct_geneIDs, by="OldPeptideID")
rownames(allresults) <- allresults$NewPeptideID
names(allresults)[names(allresults) == 'NewPeptideID'] <- 'UniqID'
names(allresults)[names(allresults) == 'New_Gene'] <- 'GeneName'
names(allresults)[names(allresults) == 'New_Pep'] <- 'UniProtID'

allresults$OldPeptideID<-NULL
allresults$Old_Gene<-NULL
allresults$Old_Pep<-NULL
allresults$t<-NULL
allresults$B<-NULL
allresults$AveExpr<-NULL

allresults$Tissue <- 'DLPFC'
names(allresults)[names(allresults) == 'logFC'] <- 'Log2_FC'
names(allresults)[names(allresults) == 'CI.L'] <- 'CI_Lwr'
names(allresults)[names(allresults) == 'CI.R'] <- 'CI_Upr'
names(allresults)[names(allresults) == 'P.Value'] <- 'PVal'
names(allresults)[names(allresults) == 'adj.P.Val'] <- 'Cor_PVal'

##reorder the columns
col_order <- c("UniqID", "GeneName", "UniProtID", "ENSG", "Tissue", "Log2_FC", "CI_Upr", "CI_Lwr", "PVal", "Cor_PVal")
allresults <- allresults[, col_order]


#save to synapse
write.csv(allresults, file="~/prot-lineage/ROSMAP_DiffExp_TMTproteins.csv", row.names=FALSE)
#save to synapse
file <- synapser::File(path='~/prot-lineage/ROSMAP_DiffExp_TMTproteins.csv', parentId='syn35219190')
file <- synapser::synStore(file)







