# RENYI ANALYSIS #

# LIBRARIES
library (BiodiversityR)
library(tidyr)
library(ggplot2)
library(extrafont) #font_import() loadfonts(device = "win") 

# DIRECTORIES
result.dir=c("C:/Users/victo/Google Drive/PAPER/counts_per_taxa_domain")
main.dir=c("C:/Users/victo/Google Drive/PAPER")
data.dir=c("C:/Users/victo/Desktop/Data") #metagenomic data

# PLOT SETTINGS
tema=theme(axis.text.x = element_text(color="black",size=15, angle=90, family="Times New Roman",hjust=0.95,vjust=0.2),
           axis.text.y = element_text(color="black",size=15, family="Times New Roman"),
           axis.title = element_text(color="black",size=15, family="Times New Roman"),
           legend.text = element_text(color = "black",size=15, family="Times New Roman"),
           panel.border =element_rect(color = "black", fill = NA) , #element_blank(),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "none")

orden=c("Am.s","Pe.s","Ma.s","Sf.s","At.Am.rz","At.Pe.rz","As.Ma.rz","As.Sf.rz", "Mg.Ma.rz","Mg.Sf.rz" ,
        "Or.Ma.rz" ,"Or.Sf.rz",  "At.Am.e","At.Pe.e","As.Ma.e","As.Sf.e","Mg.Ma.e","Mg.Sf.e", "Or.Ma.e" ,"Or.Sf.e") 

# PALETTE
aes4=c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20")
getPalette=colorRampPalette(aes4) 


# DATA (select)

#select
wplot=c("taxonomic","functional")[1]

if (wplot=="functional"){
#functional data
load(file=paste(result.dir,"taxon.ko.counts.list",sep = "/"))
data=taxon.tab.list[["Bacteria"]]+taxon.tab.list[["Archaea"]]
} else if (wplot=="taxonomic"){
#taxonomic data
data=read.delim(paste(main.dir,"family_absolut_counts.txt",sep = "/"),header = T)
rownames(data)=data$lineage ; data=data[,-1]
}

#metadata
full.meta=read.delim(paste(data.dir,"meta_data.txt",sep = "/"),header = T)
meta=data.frame(full.meta[,c("taxon_oid","species","site","compartment")], row.names = full.meta$Genome_name)


#result_norm=as.data.frame(t(t(data)/colSums(data))) #relative abundance
data=round(data, digits = 0)
datos=t(data)


##define random seed
initial.seed=as.integer(Sys.time())
the.seed=initial.seed %% 100000
print(the.seed)       #65750

the.seed=65750
set.seed(the.seed)

# Rarefied richness and alpha diversity 
ric=rarefy(datos, sample = min(colSums(data)))
datos_rar=rrarefy(datos, sample = min(colSums(data)))
rowSums(datos_rar)
shan=diversity(datos_rar,index = "shannon")


#generate plot table
result=data.frame(richness=ric, shannon=shan, row.names = names(ric))
plot.data=cbind(result[rownames(meta),],meta)
plot.data$name=rownames(plot.data)

#richness
a=ggplot(plot.data[,-2], aes(x=compartment, y=richness, fill=site, colour=site))+tema+
  geom_boxplot(fill="grey75", colour="black", size=1.3, outlier.color = "white")+
  scale_colour_manual(values = getPalette(5))+scale_fill_manual(values = getPalette(5))+
  geom_jitter(size=3, aes(shape=species), alpha=0.9)+ylab("Richness")+
  scale_x_discrete(limits=c("soil","rhizosphere","phyllosphere"),
                   labels=c("Soil","Rhizosphere","Phyllosphere"),
                   name="Compartment")+
  scale_shape_manual(values = c(21:25),
                     limits=c("No.host.associated","A.tequilana","A.salmiana","M.geometrizans","O.robusta"),
                     labels=c("Soil","A. tequilana","A. salmiana","M. geometrizans","O. robusta"))

a

png(filename = paste(result.dir,paste(wplot,"richness.png",sep = "_"),sep = "/"), 
    width = 90, height = 100, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
a ; dev.off()

#shannon index
a=ggplot(plot.data[,-1], aes(x=compartment, y=shannon, fill=site, colour=site))+tema+
  geom_boxplot(fill="grey75", colour="black", size=1.3, outlier.color = "white")+
  scale_colour_manual(values = getPalette(5))+scale_fill_manual(values = getPalette(5))+
  geom_jitter(size=3, aes(shape=species), alpha=0.9)+ylab("Shannon index")+
  scale_x_discrete(limits=c("soil","rhizosphere","phyllosphere"),
                   labels=c("Soil","Rhizosphere","Phyllosphere"),
                   name="Compartment")+
  scale_shape_manual(values = c(21:25),
                     limits=c("No.host.associated","A.tequilana","A.salmiana","M.geometrizans","O.robusta"),
                     labels=c("Soil","A. tequilana","A. salmiana","M. geometrizans","O. robusta"))

a

png(filename = paste(result.dir,paste(wplot,"shannon.png",sep = "_"),sep = "/"), 
    width = 90, height = 100, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
a ; dev.off()



#ANOSIM
the.seed=65750
datos_rar=datos_rar[as.character(full.meta$Genome_name),] #order rownames by the metadata

set.seed(the.seed)
test = anosim(datos_rar, permutations = 999, grouping = full.meta$compartment, distance = "bray")
plot(test)

set.seed(the.seed)
test = anosim(datos_rar, permutations = 999, grouping = full.meta$species, distance = "bray")
plot(test)

set.seed(the.seed)
test = anosim(datos_rar, permutations = 999, grouping = full.meta$site, distance = "bray")
plot(test)
