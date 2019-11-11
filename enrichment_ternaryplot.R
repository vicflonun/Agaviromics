### CREATE TERNARY PLOTS BASED ON KO COUNTS & KO ENRICHMENT ANALYSIS, BETWEEN COMPARTMENTS AND PLANT SPECIES ###

# SETTINGS #

# WORKING DIRECTORIES
setwd("/my_path") 
data.dir=c("/my_data") #metagenomic data
main.dir=c("/my_path/main") #enrichment analysis results
result.dir=c("/my_path/result") #saved results
dir.create(result.dir, recursive = T, showWarnings = F)

# LIBRARIES
library(data.table)
library(ggplot2)
library(ggrepel)
library(Rmisc)
library(tidyr)
library(ggtern)
library(extrafont) #font_import() loadfonts(device = "win") 

# PLOT SETTINGS
tema=theme(axis.text = element_text(color="black",size=11, angle=0, family="Times New Roman"),
           axis.title = element_text(color="black",size=11, family="Times New Roman"),
           legend.text = element_text(color = "black",size=11, family="Times New Roman"),
           panel.border =element_rect(color = "black", fill = NA) , #element_blank(),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "left") #tema 

# LOAD METADATA
path.data=fread(file = paste(getwd(),"kegg_pathway_data.txt",sep="/"),sep="\t", header = T)
pro.data=read.delim(paste(getwd(),"kegg_process_pathway_data.txt",sep = "/"), quote = "", stringsAsFactors = F)
full.meta=read.delim(paste(data.dir,"meta_data.txt",sep = "/"),header = T,stringsAsFactors = F)

# GET FILE NAMES 
files=grep(33000,list.files(data.dir),value = T) #id of samples
keys=full.meta$Genome_name ; names(keys)=full.meta$taxon_oid #match id and sample names

# LOAD COUNTS
load(file=paste(getwd(),"counts_per_taxa_domain","taxon.ko.counts.list",sep = "/"))
data.full=taxon.tab.list[["Bacteria"]]+taxon.tab.list[["Archaea"]]

# ORGANIZE DESIRED PATHWAYS
phototrophy=c("Porphyrin and chlorophyll metabolism [PATH:ko00860]","Carotenoid biosynthesis [PATH:ko00906]","Photosynthesis [PATH:ko00195]")
phosp.ox=c("Oxidative phosphorylation [PATH:ko00190]")
biofilm=c("Quorum sensing [PATH:ko02024]","Biofilm formation - Pseudomonas aeruginosa [PATH:ko02025]","Biofilm formation - Escherichia coli [PATH:ko02026]", "Flagellar assembly [PATH:ko02040]","Biofilm formation - Vibrio cholerae [PATH:ko05111]","Bacterial chemotaxis [PATH:ko02030]")
carbon.met=c("Methane metabolism [PATH:ko00680]","Glyoxylate and dicarboxylate metabolism [PATH:ko00630]","Propanoate metabolism [PATH:ko00640]","Butanoate metabolism [PATH:ko00650]")
c.fixation=c("Carbon fixation pathways in prokaryotes [PATH:ko00720]","Carbon fixation in photosynthetic organisms [PATH:ko00710]")
transporters=c("ABC transporters [PATH:ko02010]","Phosphotransferase system (PTS) [PATH:ko02060]")
two.component=c("Two-component system [PATH:ko02020]")
nitrogen.met=c("Nitrogen metabolism [PATH:ko00910]")
aminoacid=pro.data[pro.data$process=="Amino acid metabolism","pathway"]
xenobiotics=pro.data[pro.data$process=="Xenobiotics biodegradation and metabolism","pathway"]
sugar.met=c("Fructose and mannose metabolism [PATH:ko00051]","Starch and sucrose metabolism [PATH:ko00500]")
natural.products=unique(path.data[c(5518:7032),pathway])[-5] 
sulfur.met=c("Sulfur metabolism [PATH:ko00920]")

#EXTRACT THE KO IDS OF DESIRED PATHWAYS
process=list(
  Phototrophy=path.data[path.data$pathway %in% phototrophy | path.data$anygene_name %in% c("pufM","pufL","pufA","pufB","pufX","puhA","pucA","pucB"),ko_id],
  Oxidative_phosphorilation=path.data[path.data$pathway %in% phosp.ox,ko_id],
  Biofilm_QS=path.data[path.data$pathway %in% biofilm,ko_id],
  Natural_products=path.data[path.data$pathway %in% natural.products,ko_id],
  Xenobiotic_metabolism=path.data[path.data$pathway %in% xenobiotics & !(path.data$anygene_name %in% c("nifD", "nifK","nifH","anfG","vnfD","vnfK", "vnfG","vnfH")),ko_id ],
  Nitrogen_fixation=path.data[path.data$anygene_name %in% c("nifD", "nifK","nifH","anfG","vnfD","vnfK", "vnfG","vnfH") , ko_id],
  Carbon_metabolism=path.data[path.data$pathway %in% carbon.met,ko_id],
  Sugar_metabolism=path.data[path.data$pathway %in% sugar.met,ko_id],
  Carbon_fixation=path.data[path.data$pathway %in% c.fixation,ko_id],
  Nitrogen_Aminoacids=path.data[path.data$pathway %in% nitrogen.met | path.data$pathway %in% aminoacid,ko_id],
  PTS_ABC_transportes=path.data[path.data$pathway %in% transporters,ko_id],
  Two_component=path.data[path.data$pathway %in% two.component & !(path.data$anygene_name %in% c("pufM","pufL","pufA","pufB","pufX","puhA","pucA","pucB")) , ko_id]
  #sulfur.metabolism=path.data[path.data$pathway %in% sulfur.met,"ko_id"]
) #process list

#SELECT YOUR COLOR PALETTE
my_aes3=c("#024AAC","#168582","#2DBADA","#358404","#F48F06","#F2252C","#E2E41C","#10B62F","#95C331","#A30AB0", "#E283C4","hotpink")
getPalette=colorRampPalette(my_aes3) 
aes3=data.frame(color=getPalette(length(names(process))), ruta = names(process)) #color assignment


for (h in 1:2){

# PLOT BETWEEN COMPARTMENTS #

# TRANSFORM COUNTS 
result_norm=as.data.frame(t(t(data.full)/colSums(data.full))) #relative abundance
result_norm$ko_id=rownames(result_norm) #ko names
result.tidy=gather(result_norm, key = sample, value = rel.abundance, keys[files]) #create tidy data

for (i in keys){
  result.tidy[result.tidy$sample==i,c(4,5,6)]=full.meta[full.meta$Genome_name==i,c("species","site","compartment")]
} #add metadata

# SELECT DATA
file=c("Expanded_families_edgeR.txt","Reduced_families_edgeR.txt","Comparative_Genomics_families_edgeR.txt")
comparison=c("SYM_Prokaryotes_ko_75", "TEQ_Prokaryotes_ko_0.75")[h] #select the comparison 

if (sum(grep("SYM",comparison))>0){
  pair=c("phyllo_vs_rhizo","phyllo_vs_soil", "rhizo_vs_soil")
  datos=result.tidy[result.tidy$species != "A.tequilana" & !(result.tidy$site %in% c("Penjamo","Amatitan")),]
} else if (sum(grep("TEQ",comparison)>0)){
  pair=c("phyllo_vs_rhizo","phyllo_vs_soil", "rhizo_vs_soil")
  datos=result.tidy[result.tidy$species == "A.tequilana" | result.tidy$site %in% c("Penjamo","Amatitan"),]
  } #select tha data from the tidy table based on the comparison

sub1.1=read.delim(paste(main.dir,comparison,pair[1],file[1], sep = "/"))
sub1.2=read.delim(paste(main.dir,comparison,pair[1],file[2], sep = "/"))
sub2.1=read.delim(paste(main.dir,comparison,pair[2],file[1], sep = "/"))
sub2.2=read.delim(paste(main.dir,comparison,pair[2],file[2], sep = "/"))
sub3.1=read.delim(paste(main.dir,comparison,pair[3],file[1], sep = "/"))
sub3.2=read.delim(paste(main.dir,comparison,pair[3],file[2], sep = "/")) #read the expanded[1] and reduced[2] KOs per comparison [1,2,3]

# TRANSFORM DATA
mean.data=summarySE(data=datos, measurevar = "rel.abundance", groupvars = c("ko_id","compartment")) #calculate the mean relative abundances between compartmens
mean.data2=mean.data[,c(1,2,4)] #select necessary data
plot.data=spread(mean.data2, key = compartment, value = rel.abundance) #transform data from tidy to untidy in order to have the 3 axis for the triplot

# SUBSET ONLY COMPARED KOs 
all=read.delim(paste(main.dir,comparison,pair[1],file[3], sep = "/")) #read all the compared KOs
plot.data=plot.data[plot.data$ko_id %in% rownames(all),] #select them from data

# LABEL SIGNIFICANT DATA
plot.data$value="No significant"
plot.data[plot.data$ko_id %in% c(rownames(sub1.1),rownames(sub1.2),rownames(sub2.1),rownames(sub2.2),rownames(sub3.1),rownames(sub3.2)),"value"]="Significant"

# HIGHLIGHT PROCESSES
sub.process=process[names(process)[c(1:6)]]
plot.data$ruta="Other"
for (i in rev(names(sub.process))){
  plot.data[plot.data$ko_id %in% sub.process[[i]], "ruta"]=i
  } #label each KO with a process
plot.data$ruta=factor(plot.data$ruta, levels = c("Other",names(sub.process)))

# ADD FDR
all=all[plot.data$ko_id,]
plot.data$p.value=all$FDR #non significant FDR will be based on the first comparison pair[1]
p=rbind(sub1.1,sub1.2,sub2.1,sub2.2,sub3.1,sub3.2)
for (i in rownames(p)){
  plot.data[plot.data$ko_id==i,"p.value"]=p[i,"FDR"] 
  } #subtitute the FDR value of significan data of each comparison

# CREATE CHEAT DATSET TO HIGHLIGHT SPECIFIC POINTS
cheat=plot.data[plot.data$ruta %in% names(sub.process) & plot.data$value == "Significant",]

# PLOT
a=ggtern(data=plot.data, aes(x=phyllosphere, y=rhizosphere, z=soil,shape=value, color=ruta, size=-log10(p.value)))+
  geom_mask()+
  geom_point()+
  geom_point(data =cheat, aes(x=phyllosphere, y=rhizosphere, z=soil,shape=value, color=ruta,  size=-log10(p.value)))+
  scale_shape_manual(values=c(21,19))+
  scale_size_continuous(range = c(1,2))+
  scale_color_manual(values =c("grey70",as.character(aes3[aes3$ruta %in% names(sub.process),1])))+
  ylab("")+xlab("")+zlab("")+tema

#SAVE
png(filename =  paste(result.dir,paste(comparison,"ternary_9","png",sep = "."),sep = "/"), 
    width = 160, height = 100, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
print(a)  ; dev.off()


# PLOT BETWEEN SPECIES #

# TRANSFORM COUNTS
sub.meta=full.meta
sub.meta[sub.meta$species %in% c("O.robusta","M.geometrizans"),"species"]="Cacti" #merge cacti samples

result_norm=as.data.frame(t(t(data.full)/colSums(data.full))) #relative abundance 
result_norm$ko_id=rownames(result_norm) #ko names
result.tidy=gather(result_norm, key = sample, value = rel.abundance, keys[files]) #create tidy data

for (i in keys){
  result.tidy[result.tidy$sample==i,c(4,5,6)]=sub.meta[sub.meta$Genome_name==i,c("species","site","compartment")]
} #add metadata

# SELECT DATA
file=c("Expanded_families_edgeR.txt","Reduced_families_edgeR.txt","Comparative_Genomics_families_edgeR.txt")
comparison=c("RHI_Prokaryotes_ko_0.75", "PHY_Prokaryotes_ko_0.75")[h] #select the comparison

if (sum(grep("RHI",comparison))>0){
  pair=c("salmiana_vs_tequilana","salmiana_vs_cacti", "cacti_vs_tequilana")
  datos=result.tidy[result.tidy$compartment == "rhizosphere",]
} else if (sum(grep("PHY",comparison)>0)){
  pair=c("salmiana_vs_tequilana","salmiana_vs_cacti", "cacti_vs_tequilana")
  datos=result.tidy[result.tidy$compartment == "phyllosphere",]
} #select tha data from the tidy table based on the comparison

sub1.1=read.delim(paste(main.dir,comparison,pair[1],file[1], sep = "/"))
sub1.2=read.delim(paste(main.dir,comparison,pair[1],file[2], sep = "/"))
sub2.1=read.delim(paste(main.dir,comparison,pair[2],file[1], sep = "/"))
sub2.2=read.delim(paste(main.dir,comparison,pair[2],file[2], sep = "/"))
sub3.1=read.delim(paste(main.dir,comparison,pair[3],file[1], sep = "/"))
sub3.2=read.delim(paste(main.dir,comparison,pair[3],file[2], sep = "/")) #read the expanded[1] and reduced[2] KOs per comparison [1,2,3]

# TRANSFORM DATA
mean.data=summarySE(data=datos, measurevar = "rel.abundance", groupvars = c("ko_id","species")) #calculate mean relative
mean.data2=mean.data[,c(1,2,4)] #select necessary data
plot.data=spread(mean.data2, key = species, value = rel.abundance) #transform data from tidy to untidy in order to have the 3 axis for the triplot

# SUBSET ONLY COMPARED KOs 
all=read.delim(paste(main.dir,comparison,pair[1],file[3], sep = "/")) #read all the compared KOs
plot.data=plot.data[plot.data$ko_id %in% rownames(all),] #select them from data

# LABEL SIGNIFICANT DATA
plot.data$value="No significant"
plot.data[plot.data$ko_id %in% c(rownames(sub1.1),rownames(sub1.2),rownames(sub2.1),rownames(sub2.2),rownames(sub3.1),rownames(sub3.2)),"value"]="Significant"

# HIGHLIGHT PROCESSES
sub.process=process[names(process)[c(1:6)]]
sub.aes=aes3[aes3$ruta%in%names(sub.process),]

plot.data$ruta="Other"
for (i in rev(names(sub.process))){
  plot.data[plot.data$ko_id %in% sub.process[[i]], "ruta"]=i
} # label each KO with a process
plot.data$ruta=factor(plot.data$ruta, levels = c("Other",names(sub.process)))

# ADD FDR 
all=all[plot.data$ko_id,]
plot.data$p.value=all$FDR #non significant FDR will be based on the first comparison (pair[1])
p=rbind(sub1.1,sub1.2,sub2.1,sub2.2,sub3.1,sub3.2)
for (i in rownames(p)){
  plot.data[plot.data$ko_id==i,"p.value"]=p[i,"FDR"] 
} #subtitute the FDR value of significan data of each comparison

# CREATE CHEAT DATSET TO HIGHLIGHT SPECIFIC POINTS
cheat=plot.data[plot.data$ruta %in% names(sub.process) & plot.data$value == "Significant",]
cheat2=plot.data[plot.data$ruta %in% names(sub.process)[6] & plot.data$value == "Significant",]

# PLOT
a=ggtern(data=plot.data,aes(x=A.salmiana, y=A.tequilana, z=Cacti, shape=value, color=ruta, size=-log10(p.value)))+
  geom_mask()+
  geom_point()+
  geom_point(data =cheat, aes(x=A.salmiana, y=A.tequilana, z=Cacti, shape=value, color=ruta,  size=-log10(p.value)))+
  geom_point(data =cheat2, aes(x=A.salmiana, y=A.tequilana, z=Cacti, shape=value, color=ruta,  size=-log10(p.value)))+
  scale_shape_manual(values=c(21,19))+
  scale_size_continuous(range = c(1,2))+
  scale_color_manual(values =c("grey70",as.character(aes3[aes3$ruta %in% names(sub.process),1])))+
  ylab("")+xlab("")+zlab("")+tema

# SAVE
png(filename =  paste(result.dir,paste(comparison,"ternary_9","png",sep = "."),sep = "/"), 
    width = 160, height = 100, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
print(a) ; dev.off()

}

