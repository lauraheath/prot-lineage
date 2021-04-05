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
Normalized <- TMT_Express_Load('syn21266453', 1)
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


## Winzorize Expression Data
# sink_preWinzor<- Log2_Normalized
# for( i in 1:dim(Log2_Normalized)[1] ){
#   Log2_Normalized[i,] <- DescTools::Winsorize( as.numeric(Log2_Normalized[i,]), na.rm = TRUE ) 
# }
row.names(Meta) <- Meta$batchChannel


max(Log2_Normalized, na.rm=TRUE)
max(sink_preWinzor, na.rm=TRUE)
min(Log2_Normalized, na.rm=TRUE)
min(sink_preWinzor, na.rm=TRUE)

m1 <- data.matrix(Log2_Normalized)
hist(m1)









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
