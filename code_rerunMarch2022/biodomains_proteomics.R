# This script runs GSEA enrichment using branch-specific pseudotime effect sizes to rank genes and plots the results


# setup -------------------------------------------------------------------


# libraries
library(synapser)
library(fgsea)
library(tidyverse)


# plotting theme
theme_set(theme_bw())


#This function quantifies the number of significantly enriched GO terms by AD biological domain
bd.tally <- function( enrVct, biodomDefTbl){
  bdt <- bind_cols(
    domain = unique(biodomDefTbl$Biodomain),
    n_term = map_dbl( unique(biodomDefTbl$Biodomain),
                      ~ biodomDefTbl %>% filter(Biodomain == .x) %>% 
                        dplyr::select(GOterm_Name) %>% distinct() %>% nrow()),
    n_sig_term = map_dbl( unique(biodomDefTbl$Biodomain),
                          ~ enrVct[ enrVct %in% 
                                      biodomDefTbl$GOterm_Name[
                                        biodomDefTbl$Biodomain == .x]] %>% 
                            length()) ) %>% 
    mutate(domain = fct_reorder(domain, n_sig_term, .desc = F)) %>% 
    arrange(domain)
  bdt$prop <- bdt$n_sig_term / bdt$n_term
  return(bdt)
}


# data: biodomain definitions, biodomain plotting colors, and braak state-gene table
biodom <- readRDS( synapser::synGet('syn25428992')$path  )
biodom.annotated <- biodom %>% filter(!is.na(n_symbol)) %>% pull(symbol, name=GOterm_Name)
dom.cols <- read_csv( synGet('syn26856828')$path )

#pt_data <- read_csv( synGet('syn27260460')$path )
#pt_data <- read.csv(file="~/neuropath_lineages/proteomics/Dx_analysis/female_DEanova_stats.csv")
pt_data <- read.csv(file="~/neuropath_lineages/proteomics/Dx_analysis/male_DEanova_stats.csv")

# Run enrichments by pseudotime state -------------------------------------


# enrichment for each state; all significant genes (pval < 0.05), arranged by effect size
enr <- map_dfr(
  unique(pt_data$state),
  ~ pt_data %>% 
    filter(state == .x
           , pvalue < 0.05
           # , !(gene_names %in% (dlpfc.pt %>% filter(pvalue <= 0.05) %>% pull(gene_names))) # use this to consider genes unique to the neuropath
    ) %>% #
    arrange(desc(effect)) %>% 
    pull(effect, name = gene_short_name) %>%
    fgseaMultilevel(biodom.annotated, ., eps = 0, scoreType = 'std') %>% 
    mutate(state = .x)
)


# count up the number of significant terms per biodomain
bdt <- bd.tally(enr$pathway[enr$pval <= 0.05], biodom) %>% 
  mutate(domain = fct_reorder(domain, n_sig_term, .desc = T)) %>% 
  arrange(domain)


# add biodomain annotations and plotting colors to enrichment results
enr <- enr %>% 
  left_join(., biodom %>% select(pathway=GOterm_Name, Biodomain), by='pathway') %>%
  full_join(., dom.cols, by = c('Biodomain'='domain')) %>%
  mutate( Biodomain = fct_relevel(Biodomain, as.character(bdt$domain) )) %>% 
  mutate( state = fct_reorder( as.factor(state), as.numeric(state), .desc = T ))


# plot! -------------------------------------------------------------------


enr %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enr$pval) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enr, pval > 0.05),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enr, pval < 0.05 ),
    aes(color = color, size = -log10(pval) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enr, pval < 0.05 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enr, pval < 0.05 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Biodomain)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  ggtitle('Male proteomics, AD case/control DE genes; state pval < 0.05; term pval < 0.05')


enr$leadingEdge <- as.character(enr$leadingEdge)

write.csv(enr, file="~/neuropath_lineages/proteomics/Dx_analysis/F_GOenrichment.csv", row.names=FALSE)
file <- synapser::File(path='~/neuropath_lineages/proteomics/Dx_analysis/F_GOenrichment.csv', parentId='syn28558769')
file <- synapser::synStore(file)


write.csv(enr, file="~/neuropath_lineages/proteomics/Dx_analysis/M_GOenrichment.csv", row.names=FALSE)
file <- synapser::File(path='~/neuropath_lineages/proteomics/Dx_analysis/M_GOenrichment.csv', parentId='syn28558769')
file <- synapser::synStore(file)
