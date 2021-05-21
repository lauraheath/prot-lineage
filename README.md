# prot-lineage
code to perform lineage analysis on proteomic data

Requirements: users need to be a an AMP-AD Consortium member on Synapse to pull all data. Synapse IDs required are embedded in the code.
This code will use the Monocle2 R package to calculate pseudotimes based on TMT-proteomic changes from ROSMAP DLPFC samples (N=400).
All analyses are sex-stratified, see comments in code.

To obtain pseudotimes and differentially expressed genes by branch, run following R scripts in order:
1. packages_dependencies.r
2. get_proteins_metadata.r
3. protein_lineage_monocle2.R
4. DE_GOenrichment.R

To group state-specific GO enrichment terms by biodomain (annotated by Jesse Wiley, R markdown notebook by Greg Cary):
5. pseudotimes_biodoms.Rmd

To make Venn diagrams and see correlation between RNA-seq pseudotimes and proteomics pseudotimes:
6. rerun_rnaseq_lineage.r
7. prot_vs_bulk_comparisons.r


 
