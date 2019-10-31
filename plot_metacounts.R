### PLOT THE KO CATEGORIES IN EACH METAGENOME BY LINEAGE AND PERFORMS A NMDS##

# WORKING DIRECTORIES
setwd("C:/my_path")
data.dir=c("C:/my_data") #metagenomic data
result.dir=c("C:/my_path/counts_per_taxa_domain")
dir.create(result.dir, recursive = T, showWarnings = F)

# LIBRARIES
library(data.table)
library(Rmisc)
library(splitstackshape)
library(tidyr)
library(ggplot2)
library(extrafont) #font_import() loadfonts(device = "win") 

# PLOT SETTINGS
tema=theme(axis.text.x = element_text(color="black",size=13, angle=90, family="Times New Roman",hjust=0.95,vjust=0.2),
           axis.text.y = element_text(color="black",size=13, family="Times New Roman"),
           axis.title = element_text(color="black",size=13, family="Times New Roman"),
           legend.text = element_text(color = "black",size=13, family="Times New Roman"),
           panel.border =element_blank(), #element_rect(color = "black", fill = NA) ,
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "none")

orden=c("Am.s","Pe.s","Ma.s","Sf.s","At.Am.rz","At.Pe.rz","As.Ma.rz","As.Sf.rz", "Mg.Ma.rz","Mg.Sf.rz" ,
        "Or.Ma.rz" ,"Or.Sf.rz",  "At.Am.e","At.Pe.e","As.Ma.e","As.Sf.e","Mg.Ma.e","Mg.Sf.e", "Or.Ma.e" ,"Or.Sf.e") 

# PALETTE
aes4=c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20")
getPalette=colorRampPalette(aes4) 

# SET LEVELS AND COLUMN NAMES OF DATA
taxa.col=c("gene_id","homolog_gene_oid", "homolog_taxon_oid", "percent_identity","taxonomy")
depth.col=c("orig.id","depth_value")

# GET FILE NAMES 
files=grep(33000,list.files(data.dir),value = T) #id of samples
full.meta=read.delim(paste(data.dir,"meta_data.txt",sep = "/"),header = T,stringsAsFactors = F) #metadata
keys=full.meta$Genome_name ; names(keys)=full.meta$taxon_oid #match sample name and id

# LOAD PATHWAY COUNTS AND METADATA
process=read.delim(file=paste(getwd(),"kegg_process_pathway_data.txt",sep="/"))
process=process[order(process$pathway),] #order process by the order of pathways
load(file=paste(result.dir,"taxon.paths.counts.list",sep = "/"))

# SELECT THE TAXA TO WORK
data=p.taxon.tab.list[["Bacteria"]]+p.taxon.tab.list[["Archaea"]]
#i=taxon[2]
#data=p.taxon.tab.list[[i]]

# TRANSFORM DATA
result_norm=as.data.frame(t(t(data)/colSums(data))) #relative abundance
result_norm$path=rownames(result_norm) #create column with path names
result_norm$process=process$process #create column with processes
result.tidy=gather(data=result_norm, key = assembled.library, value = relative.abundance, keys[files]) #tidy format

# COLORS
colourCount = length(unique(result.tidy$process))

# PLOT
a=ggplot(data=result.tidy, aes(x=assembled.library, y=relative.abundance*100, fill=process))+
  geom_bar(stat = "identity",width = 0.8)+
  scale_fill_manual(values = getPalette(colourCount))+tema+
  scale_x_discrete(limits=orden, labels=orden)+ylab("Relative abundance (%)")

a

# SAVE plot
png(filename = paste(result.dir,"relative_abundance_processes.png",sep = "/"), 
    width = 130, height = 100, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
a ; dev.off()

#save legend
png(filename = paste(result.dir,"relative_abundance_processes_legend.png",sep = "/"), 
    width = 350, height = 160, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
a ; dev.off()



# NMDS OF KO COUNTS#

library(vegan)
library(labdsv)
library(MASS)
library(ggrepel)
library(RColorBrewer)

# LOAD DATA
load(file=paste(result.dir,"taxon.ko.counts.list",sep = "/"))
data=taxon.tab.list[["Bacteria"]]+taxon.tab.list[["Archaea"]]
result_norm=as.data.frame(t(t(data)/colSums(data))) #calculate relative abundance

result_norm=result_norm[,-9] #eliminate very dismilar samples

# CALCULATE DIMENSIONS
scaling=vegdist(t(result_norm), method = "bray") #calculate distances
scaling2=isoMDS(scaling) #calculate NMDS
#plot(hclust(scaling, method = "average")) #plot dendogram
scaling3=scaling2$points #get points

scaling3=as.data.frame(scaling3[full.meta$Genome_name[full.meta$Genome_name %in% colnames(result_norm)],])
scaling3$plant=full.meta$species[full.meta$Genome_name %in% colnames(result_norm)]
scaling3$comp=full.meta$compartment[full.meta$Genome_name %in% colnames(result_norm)]
scaling3$sample=rownames(scaling3) #add metadata
scaling3$comp=factor(scaling3$comp, levels = c("soil","rhizosphere","phyllosphere")) #order levels

# PLOT
a=ggplot(data=scaling3, aes(x=scaling3$V1, y=scaling3$V2, fill=comp, labels=sample, shape=plant))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=4, alpha=0.9)+
  #geom_text_repel(label=scaling3$sample, size=3, family="Times New Roman", color="black")+
  scale_fill_manual(values = getPalette(20)[c(1,13,5)])+tema+
  scale_colour_manual(values = getPalette(20)[c(1,13,5)])+
  scale_shape_manual(values= c(21:25),
                     limits=c("No.host.associated","A.tequilana","A.salmiana","M.geometrizans","O.robusta"),
                     labels=c("Soil","A. tequilana","A. salmiana","M. geometrizans","O. robusta"))+
  xlab("Axis_1")+ylab("Axis_2")+scale_y_continuous(limits = c(-0.7,0.7),
                                                   breaks = c(-0.7,-0.3,0,0.3,0.7))

a

# SAVE
png(filename =  paste(result.dir,"KO.nmds_final.png",sep = "/"), 
    width = 80, height = 80, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
a ; dev.off()

