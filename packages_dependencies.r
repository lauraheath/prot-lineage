
#UPDATED:
#dependencies to run in terminal window before installing packages:
#sudo apt-get update -y
#sudo apt-get install -y libudunits2-dev libgdal-dev

## also set up synapse login configuration in terminal:
#nano .synapseConfig
##in nano window type the following:
#[authentication]
#username = <username>
#password = <password>

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

#install old version of RcppAnnoy in order to successfully install later packages from bioconductor:
packageurl <- "https://cran.r-project.org/src/contrib/Archive/RcppAnnoy/RcppAnnoy_0.0.14.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

BiocManager::install(c("Biobase", "BiocGenerics","DelayedArray","DelayedMatrixStats","limma","S4Vectors","SingleCellExperiment","SummarizedExperiment", "batchelor", "Matrix.utils"))
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
BiocManager::install("monocle")
BiocManager::install("biomaRt")

install.packages("kableExtra")
install.packages("DescTools")
install.packages("ggpubr")
install.packages("enrichR")

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(monocle3)
library(monocle)
library(limma)
library(Matrix)
library(kableExtra)
library(DescTools)
library(biomaRt)
library(enrichR)

