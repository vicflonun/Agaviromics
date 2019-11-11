### CREATES BARPLOTS AND NON-METRIC MULTIDIMENSIONAL SCALING (NMDS) OF GENECOUNTS BASED ON TAXONOMIC ANNOTATION ###

# WORKING DIRECTORIES
setwd("C:/my_path")
data.dir=c("C:/my_data") #metagenomic data
main.dir=c("C:/my_path/main/")
result.dir=c("C:/my_path/main/result")
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

# PLOT TAXONOMIC RANKS

#data
taxd=c("domain","phylum","class","order","family","genus")[3] #select rank
load(file=paste(main.dir,paste("meta",taxd,"ids",sep = "_"),sep="/")) #load all phylum names
result=read.delim(file = paste(main.dir,paste(taxd,"absolut_counts.txt",sep="_"),sep="/"),sep="\t", header = T) #load data
result_norm=as.data.frame(t(t(result[,-1])/colSums(result[,-1]))*100) #relative abundance
result_norm$lineage=result$lineage #replace lineage names
result.tidy=gather(data=result_norm, key = assembled.library, value = relative.abundance, keys[files]) #transform to tidy format

#order
lin.orden=summarySE(result.tidy, measurevar = "relative.abundance", groupvars = c("lineage")) #mean abundance in all metagenomes
lin.orden=as.character(lin.orden[order(lin.orden$relative.abundance,decreasing = T),"lineage"]) #sort

#filter by abundance
filter=1
result.tidy$lineage=factor(result.tidy$lineage, levels = c(lin.orden,paste("Lineages<",filter,"%",sep = ""))) #set levels with the filter
result.tidy[result.tidy$relative.abundance<filter,"lineage"]=paste("Lineages<",filter,"%",sep = "") #filter
ncolor = length(unique(result.tidy$lineage));print(ncolor) #sanity check for the numer of taxa (better less than 12)

#plot
a=ggplot(data=result.tidy, aes(x=assembled.library, y=relative.abundance, fill=lineage))+
  geom_bar(stat = "identity",width = 0.8)+
  scale_fill_manual(values = getPalette(ncolor))+tema+
  scale_x_discrete(limits=orden, labels=orden)+ylab("Relative abundance (%)")+xlab("Sample")+
  guides(col = guide_legend(ncol= 2))

#save plot
png(filename = paste(result.dir,paste(taxd,"relative_abundance_plot.png",sep="_"),sep = "/"), 
    width = 130, height = 100, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
a ; dev.off()

#save legend
png(filename = paste(result.dir,paste(taxd,"relative_abundance_legend.png",sep="_"),sep = "/"), 
    width = 350, height = 160, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
a ; dev.off()



# NMDS #

library(vegan)
library(labdsv)
library(MASS)
library(ggrepel)
library(RColorBrewer)

#data
taxd=c("domain", "phylum", "class","order","family","genus")[4]
result=read.delim(file = paste(main.dir,paste(taxd,"absolut_counts.txt",sep="_"),sep="/"),sep="\t", header = T,stringsAsFactors = F)
result_norm=as.data.frame(t(t(result[,-1])/colSums(result[,-1])))
rownames(result_norm)=result$lineage

#calculate dimensions
scaling=vegdist(t(result_norm), method = "bray") #calculate distance
scaling2=isoMDS(scaling) ; scaling2$stress #create NMDS 
scaling3=scaling2$points #select cordenates
scaling3=as.data.frame(scaling3[full.meta$Genome_name,]) #aort as metadata
scaling3$plant=full.meta$species
scaling3$comp=full.meta$compartment
scaling3$sample=rownames(scaling3) #add metadata
scaling3$comp=factor(scaling3$comp, levels = c("soil","rhizosphere","phyllosphere")) #define order

#plot
a=ggplot(data=scaling3, aes(x=scaling3$V1, y=scaling3$V2, fill=comp, labels=sample, shape=plant))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=4, alpha=0.9)+
  #geom_text_repel(label=scaling3$sample, size=3, family="Times New Roman", color="black")+
  #scale_colour_manual(values = getPalette(20)[c(1,13,5)])+
  scale_shape_manual(values= c(21:25),
                     limits=c("No.host.associated","A.tequilana","A.salmiana","M.geometrizans","O.robusta"),
                     labels=c("Soil","A. tequilana","A. salmiana","M. geometrizans","O. robusta"))+
  scale_fill_manual(values = getPalette(20)[c(1,13,5)])+tema+
  xlab("Axis_1")+ylab("Axis_2")+scale_y_continuous(limits = c(-0.7,0.7),
                                                   breaks = c(-0.7,-0.3,0,0.3,0.7))

a

#save
png(filename =  paste(result.dir,paste(taxd,"nmds_final.png", sep = "_"),sep = "/"), 
    width =80, height = 80, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
a ; dev.off()

