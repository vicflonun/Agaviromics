# Agaviromics
This repository contains the code used to perform the downstream statistical analyses in Flores-Nuñez, et al. submmited to Frontiers in Microbiology. 

get_taxacounts.R : Generates the taxonomic profile of the genes of each library, from domain to genus
get_metacounts.R : Generates the functional profile of the genes of each library based on KO annotation 
meta_diversity.R : Calculates diversity indexes from raified taxonomic and functional counts
plot_taxacounts.R : Creates barplots and NMDS the taxonomic gene counts
plot_metacounts.R : Creates barplots and NMDS the functional (KO) gene counts
meta_enrichment.R : Performs the enrichment analysis between different subset of grouped samples using edgeR
meta_ova.R : Performs an overrepresentation anaylis of the pathways (KEGG) in each gene enrichment list
meta_AAP.R : Creates abundance plots of aerobic anoxygenic phototrophy related genes
meta_nif.R : Creates abundance plots of nif genes
meta_biofilm.R : Creates abundance plots of biofilm and QS genes
