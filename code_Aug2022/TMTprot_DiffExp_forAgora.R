library(limma)

## Import batch-corrected, centered, Log2-transformed protein expression matrix
p <- synapser::synGet('syn21266454')
Log2_Normalized <- read.csv(p$path)

#Fix protein names (some out of date)
p1 <- synapser::synGet('syn24216770')
correct_geneIDs <- read.csv(p1$path)

names(Log2_Normalized)[names(Log2_Normalized) == 'X'] <- 'OldPeptideID'
Log2_Normalized <- dplyr::left_join(Log2_Normalized, correct_geneIDs, by="OldPeptideID")
rownames(Log2_Normalized) <- Log2_Normalized$NewPeptideID
Log2_Normalized$OldPeptideID<-NULL
Log2_Normalized$NewPeptideID<-NULL
Log2_Normalized$Old_Gene<-NULL
Log2_Normalized$Old_Pep<-NULL
Log2_Normalized$New_Gene<-NULL
Log2_Normalized$New_Pep<-NULL
Log2_Normalized$ENSG<-NULL

#get sample metadata
p2 <- synapser::synGet('syn21323404')
Meta <- read.csv(p2$path)

# - get patient BioSpecimen Data with rosmap individual ID
p3 <- synapser::synGet('syn21323366')
BioSpecimen <- read.csv(p3$path)
BioSpecimen <- BioSpecimen[ BioSpecimen$assay == 'TMT quantitation', ]
#get rid of pooled samples
Meta <- Meta[ Meta$isAssayControl==FALSE,]
Meta <- dplyr::left_join(Meta, BioSpecimen, by = 'specimenID' )

#get clinical metadata
p4 <- synapser::synGet('syn3191087')
Clinical <- read.csv(p4$path)
Meta <- dplyr::left_join(Meta, Clinical, by = 'individualID' )
Meta <- Meta[ ,colSums(is.na(Meta))<nrow(Meta) ] 
row.names( Meta ) <- Meta$batchChannel
Meta <- Meta[ colnames(Log2_Normalized), ]

Meta$braaksc <- as.numeric(Meta$braaksc)
Meta$ceradsc <- as.numeric(Meta$ceradsc)
Meta$cogdx <- as.numeric(Meta$cogdx)
# Harmonize case-control status
Meta$diagnosis <- "other"
Meta$diagnosis[Meta$cogdx == 1 & Meta$braaksc <= 3 & Meta$ceradsc >= 3] <- "control"
Meta$diagnosis[Meta$cogdx == 4 & Meta$braaksc >= 4 & Meta$ceradsc <= 2] <- "AD"
table(Meta$diagnosis)

## Add Ages over 90 for modeling
p5 <- synapser::synGet('syn23573928')
Mast <- read.csv(p5$path)
Mast<- Mast[ !duplicated(Mast$individualID), ]
Meta <- dplyr::left_join(Meta[ , colnames(Meta)[ (colnames(Meta) %in% c('age_at_visit_max', 'age_first_ad_dx', 'age_death') )==F ] ], Mast[, c('individualID', 'age_at_visit_max', 'age_first_ad_dx', 'age_death')],by = 'individualID')
#Convert APOE into numeric (0,1,2 for #e4 alleles)
Meta$apoe_genotype <- as.numeric( Meta$apoe_genotype )
APOE <- as.numeric(names( table(Meta$apoe_genotype) ))
names(APOE) <- APOE
APOE[names( table(Meta$apoe_genotype) )] <- 0
APOE[grepl("4",names(APOE))] <- 1
APOE[grepl("44",names(APOE))] <- 2
Meta$APOE <- as.factor( APOE[ as.character(Meta$apoe_genotype) ] )

####Metadata file for AD case-control Diagnosis only for analysis
Meta_D <- Meta[ Meta$diagnosis %in% c('AD','control'), ]
#Code Sex and Diagnosis as Factors
Meta$msex <- as.factor(Meta$msex)
Meta_D$diagnosis <- as.factor(Meta_D$diagnosis)
#pare log2 matrix to AD case and controls
ADmatrix <- Log2_Normalized[Meta_D$batchChannel]

#run limma, adjust for sex and pmi 
design <- model.matrix(~ 0 + diagnosis + msex + pmi, data = Meta_D)
design
cont.matrix <- makeContrasts(diagnosisAD-diagnosiscontrol, levels=colnames(design))
cont.matrix

fit <- lmFit(ADmatrix, design)
fit <- contrasts.fit(fit, cont.matrix)
#do robust regression to minimize effect of outlier genes
fit <- eBayes(fit, robust=TRUE)
topTable(fit)

#extract results for all genes and get 95% confidence limits
allresults <- topTable(fit, n=Inf, confint=0.95)

#relabel columns to match lfq data for agora, and get ENSG names

allresults$NewPeptideID <- rownames(allresults)
allresults <- dplyr::left_join(allresults, correct_geneIDs, by="NewPeptideID")
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
write.csv(allresults, file="~/prot-lineage/data_objects/ROSMAP_DiffExp_TMTproteins.csv", row.names=FALSE)
#save to synapse
file <- synapser::File(path='~/prot-lineage/data_objects/ROSMAP_DiffExp_TMTproteins.csv', parentId='syn35219190')
file <- synapser::synStore(file)

#save data objects to local folder for downstream pseudotime analyses
write.csv(Log2_Normalized, file="~/prot-lineage/data_objects/Log2_Normalized.csv")
write.csv(Meta, file="~/prot-lineage/data_objects/TMT_metadata.csv")





