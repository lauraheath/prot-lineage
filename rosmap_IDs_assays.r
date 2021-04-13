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

Meta2 <- merge(rosmap_rnaseq, Meta, all=TRUE)
Meta2$rnaseq[is.na(Meta2$rnaseq)] <- 0
Meta2$proteomics[is.na(Meta2$proteomics)] <- 0
table(Meta2$rnaseq, Meta2$proteomics)
Meta2 <- subset(Meta2, Meta2$proteomics==1)



#Code Sex and Diagnosis as Factors
Meta$msex <- as.factor(Meta$msex)
Meta$APOE <- as.factor(Meta$APOE)

####Meta for Diagnosis
Meta_D <- Meta[ Meta$diagnosis %in% c('AD','control'), ]
Meta_D$diagnosis <- as.factor(Meta_D$diagnosis)


##### for patient characteristics table ####
Meta$braaksc <- as.factor(Meta$braaksc)
Meta$ceradsc <- as.factor(Meta$ceradsc)
Meta$cogdx <- as.factor(Meta$cogdx)

table(Meta$msex)
prop.table(table(Meta$msex))

sum(is.na(Meta$diagnosis))
table(Meta$diagnosis, Meta$msex)
female <- subset(Meta, Meta$msex==0)
prop.table(table(female$diagnosis))
male <- subset(Meta, Meta$msex==1)
prop.table(table(male$diagnosis))

mean(female$age_death, na.rm=TRUE)
sd(female$age_death, na.rm=TRUE)
mean(male$age_death, na.rm=TRUE)
sd(male$age_death, na.rm=TRUE)

mean(female$age_death, na.rm=TRUE)
sd(female$age_death, na.rm=TRUE)
mean(male$age_death, na.rm=TRUE)
sd(male$age_death, na.rm=TRUE)

table(female$braaksc)
round(prop.table(table(female$braaksc)), digits=3)
sum(is.na(female$braaksc))
table(male$braaksc)
round(prop.table(table(male$braaksc)), digits=3)
sum(is.na(male$braaksc))

table(female$ceradsc)
round(prop.table(table(female$ceradsc)), digits=3)
sum(is.na(female$ceradsc))
table(male$ceradsc)
round(prop.table(table(male$ceradsc)), digits=3)
sum(is.na(male$ceradsc))


table(female$cogdx)
round(prop.table(table(female$cogdx)), digits=3)
sum(is.na(female$cogdx))
table(male$cogdx)
round(prop.table(table(male$cogdx)), digits=3)
sum(is.na(male$cogdx))

table(female$apoe_genotype)
round(prop.table(table(female$apoe_genotype)), digits=3)
sum(is.na(female$apoe_genotype))
table(male$apoe_genotype)
round(prop.table(table(male$apoe_genotype)), digits=3)
sum(is.na(male$apoe_genotype))



table(Meta$apoe_genotype, Meta$msex)
table(Meta$ceradsc, Meta$msex)
table(Meta$cogdx, Meta$msex)



### get number of people in TMT who were also in the bulk rnaseq lineage paper, by sex

Meta3 <- subset(Meta2, Meta2$proteomics==1)
table(Meta3$rnaseq, Meta3$msex)

### jake's code with various IDs ########
synapser::synLogin()
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
#Indv MetaData
Indv <- TMT_Express_Load('syn3191087', 0)
pro_Assay <- TMT_Express_Load('syn23569441', 0)
tran_Assay <- TMT_Express_Load('syn21088596', 0)
Bio <- TMT_Express_Load('syn21323366', 0)
Meta <- dplyr::left_join(pro_Assay, Bio, by='specimenID')
Meta_b <- dplyr::left_join(tran_Assay, Bio, by='specimenID')
META <- dplyr::left_join(Meta, Indv, by='individualID')
META_final <- dplyr::left_join(Meta_b, META, by='individualID')
