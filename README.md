# Agaviromics
This repository contains the code used to perform the downstream statistical analyses in Flores-Nuñez, et al. 2020 in Frontiers in Microbiology. 

The datasets generated for this study can be found in the IMG platform, IMG Genome IDs 3300032159,
3300030499, 3300030692, 3300030502, 3300030753, 3300030500, 3300030504, 3300030512, 3300030505, 3300030514, 3300030501,
3300030497, 3300030498, 3300030516, 3300030495, 3300030515, 3300030511, 3300030510, 3300030496, 3300030513.

# get_taxacounts.R : Generates the taxonomic profile of the genes of each library, from domain to genus
# get_metacounts.R : Generates the functional profile of the genes of each library based on KO annotation 
# meta_diversity.R : Calculates diversity indexes from raified taxonomic and functional counts
# plot_taxacounts.R : Creates barplots and NMDS the taxonomic gene counts
# plot_metacounts.R : Creates barplots and NMDS the functional (KO) gene counts
# meta_enrichment.R : Performs the enrichment analysis between different subset of grouped samples using edgeR
# meta_ova.R : Performs an overrepresentation anaylis of the pathways (KEGG) in each gene enrichment list
# meta_AAP.R : Creates abundance plots of aerobic anoxygenic phototrophy related genes
# meta_nif.R : Creates abundance plots of nif genes
# meta_biofilm.R : Creates abundance plots of biofilm and QS genes
