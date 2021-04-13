synapser::synLogin()

TMT_Express_Load <- function( SynID, names ){
  #'@SynID a synapse ID of a proteomics csv matrix eg. syn21266454 or syn21266454
  #'@names number corresponding to the row that are the rownames (Use 0 if none apply )
  if( !is.character(SynID) ){
    return("Invalid Synapse ID: Input SynID is not a character string")
  }
  if( !grepl('syn', SynID) | !((nchar(SynID) == 11) | (nchar(SynID) == 10 ))  ){
    return("Invalid Synapse ID: Input SynID is not a valid synapse ID")
  }
  #if(  grepl('Benefactor not found for', as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) )) ){
  #  return("Syn ID does not exist")
  #}
  if(  grepl('You do not have READ permission for the requested entity', as.character( try( synapser::synGetPermissions(SynID), silent = TRUE ) )) ){
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

Log2_Normalized <-TMT_Express_Load('syn21266454', 1)
#Normalized <- TMT_Express_Load('syn21266453', 1)
Meta <- TMT_Express_Load('syn21323404', 1)

# - Public Facing BioSpecimin Data: syn21323366
# - Staged BioSpecimin Data: syn23583548

p2 <- synapser::synGet('syn23583548')
BioSpecimin <- read.csv(p2$path)

#BioSpecimin <- TMT_Express_Load('syn21323366', 0)

#if( "assay" %in% colnames(BioSpecimin) ){
#}else{
#  BioSpecimin <- TMT_Express_Load('syn23583548', 0)
#}

BioSpecimin <- BioSpecimin[ BioSpecimin$assay == 'TMT quantitation', ]
Meta$specimenID <- row.names(Meta)
Meta <- dplyr::left_join(Meta, BioSpecimin, by = 'specimenID' )

Clinical <- TMT_Express_Load( 'syn3191087',0 )

Meta <- dplyr::left_join(Meta, Clinical, by = 'individualID' )
Meta <- Meta[ ,colSums(is.na(Meta))<nrow(Meta) ] 
row.names( Meta ) <- Meta$batchChannel
Meta <- Meta[ colnames(Log2_Normalized), ]
Meta <- Meta[ ,colnames(Meta)[ (colnames(Meta) %in%'controlType' )==F] ]

# Harmonize case-control status
Meta$diagnosis <- "Other"
Meta$diagnosis[Meta$cogdx == 1 & Meta$braaksc <= 3 & Meta$ceradsc >= 3] <- "Control"
Meta$diagnosis[Meta$cogdx == 4 & Meta$braaksc >= 4 & Meta$ceradsc <= 2] <- "AD"

kableExtra::kable(table(Meta$diagnosis))

## Add Ages over 90 for modeling
#Mast <- TMT_Express_Load('syn23573928', 0)
p <- synapser::synGet('syn23573928')
Mast <- read.csv(p$path)
Mast<- Mast[ !duplicated(Mast$individualID), ]
Meta <- dplyr::left_join(Meta[ , colnames(Meta)[ (colnames(Meta) %in% c('age_at_visit_max', 'age_first_ad_dx', 'age_death') )==F ] ], Mast[, c('individualID', 'age_at_visit_max', 'age_first_ad_dx', 'age_death')],by = 'individualID')

#Convert APOE
APOS <- as.numeric(names( table(Meta$apoe_genotype) ))
names(APOS) <- APOS
APOS[names( table(Meta$apoe_genotype) )] <- 0
APOS[grepl("4",names(APOS))] <- 1
APOS[grepl("44",names(APOS))] <- 2
Meta$APOE <- as.numeric( APOS[ as.character(Meta$apoe_genotype) ] )

#impute missing pmi (2 missing: fill in with median pmi)
paste('Imputing PMI to:',median(Meta$pmi[!is.na(Meta$pmi)]))
#add this back into the metadata file
Meta$pmi[is.na(Meta$pmi)] <- 6.5

row.names(Meta) <- Meta$batchChannel

#add an indicator for whether the sample was also in the bulk rna-seq lineage analysis. need to find common ID:
dlpfcCovObj <- synapser::synGet('syn11024258')
rosmap1 <- read.delim(dlpfcCovObj$path,stringsAsFactors = F)
rosmapObj <- synapser::synGet('syn3191087')
rosmap2 <- data.table::fread(rosmapObj$path,data.table=F)

#removing bad batches in rosmap
rosmap1 <- subset(rosmap1, rosmap1$Batch<7)


#need to synchronize Sample IDs & join)
rosmapIdObj <- synapser::synGet('syn3382527')
rosmapId <- data.table::fread(rosmapIdObj$path,data.table=F)
rosmapId <- dplyr::select(rosmapId,projid,rnaseq_id)
rosmapRNAid<-dplyr::left_join(rosmapId,rosmap2)
#remove duplicate rows
rosmapRNAid <- unique(rosmapRNAid)

rosmapRNAid2 <- subset(rosmapRNAid, select=c(rnaseq_id, projid, individualID))
#names(rosmapRNAid2)[names(rosmapRNAid2) == "rnaseq_id"] <- "SampleID"
names(rosmap1)[names(rosmap1) == "SampleID"] <- "rnaseq_id"

rosmap_rnaseq <-dplyr::left_join(rosmap1,rosmapRNAid2, by="rnaseq_id")
rosmap_rnaseq <- subset(rosmap_rnaseq, select=c(rnaseq_id, individualID))
rosmap_rnaseq$rnaseq = 1

Meta$proteomics = 1

Meta <- merge(rosmap_rnaseq, Meta, all=TRUE)
Meta$rnaseq[is.na(Meta$rnaseq)] <- 0
Meta$proteomics[is.na(Meta$proteomics)] <- 0
table(Meta$rnaseq, Meta$proteomics)
Meta <- subset(Meta, Meta$proteomics==1)

rownames(Meta) <- Meta$batchChannel

#need to reorder Meta rows to match the protein data column order
meta_ordered <- match(colnames(Log2_Normalized), rownames(Meta))
Meta <- Meta[meta_ordered,]


#fix gene | uniprot identifiers (per Jake: some genes changed after checking ensembl gene ids)
p <- synapser::synGet('syn24216770')
correct_geneIDs <- read.csv(p$path)
Log2_Normalized$OldPeptideID <- rownames(Log2_Normalized)
Log2_Normalized <- dplyr::left_join(Log2_Normalized, correct_geneIDs, by="OldPeptideID")
rownames(Log2_Normalized3) <- Log2_Normalized$NewPeptideID
Log2_Normalized$Old_Gene<-NULL
Log2_Normalized$Old_Pep<-NULL
Log2_Normalized$OldPeptideID<-NULL
Log2_Normalized$New_Gene<-NULL
Log2_Normalized$New_Pep<-NULL
Log2_Normalized$NewPeptideID<-NULL
Log2_Normalized$ENSG<-NULL




# # change the '|' to a period in the row names of the protein matrices
# Log2_Normalized2 <- Log2_Normalized
# Log2_Normalized2$proteins <- rownames(Log2_Normalized2)
# Log2_Normalized2$proteins <- gsub("\\|.*", "", Log2_Normalized2$proteins)
# #need to make duplicated gene short names unique
# Log2_Normalized2$proteins <- make.names(Log2_Normalized2$proteins, unique=TRUE)
# rownames(Log2_Normalized2) <- Log2_Normalized2$proteins
# Log2_Normalized2$proteins<-NULL
# 
# 
# #check out some protein distributions
# log2prots <- t(Log2_Normalized)
# log2prots <- as.data.frame(log2prots)
# hist(log2prots$VAMP1_P23763)
# hist(log2prots$WDR33_Q9C0J8)
# sum(is.na(log2prots$KCTD13_Q8WZ19))
# sum(is.na(log2prots$CBX3_Q13185))
# max(log2prots$WDR33)
# 
# log2prots <- t(sink_preWinzor)
# log2prots <- as.data.frame(log2prots)
# hist(log2prots$VAMP1)
# hist(log2prots$WDR33)
# sum(is.na(log2prots$KCTD13))
# sum(is.na(log2prots$CBX3))
# max(log2prots$`WDR33|Q9C0J8`)
# summary(log2prots$`WDR33|Q9C0J8`)
# 
# 
# 
# rando <- sample_n(Log2_Normalized2, 20)
# rando2 <- data.matrix(rando)
# hist(Normalized)
# hist.data.frame(rando)
# 
# datmatrix <- data.matrix(Log2_Normalized)
# hist(datmatrix)
# install.packages("Hmisc")
# library(Hmisc)
# hist.data.frame(rando)
# h <- ggplot(data = log2prots, aes(x = ))
