### CREATE A TAXONOMIC PROFILE OF THE METAGENOMES ###

# WORKING DIRECTORIES
setwd("C:/my_path")
data.dir=c("C:/my_data") #metagenomic data
result.dir=c("C:/my_path/PAPER/")
dir.create(result.dir, recursive = T, showWarnings = F)

# LIBRARIES
library(data.table)
library(Rmisc)
library(splitstackshape)
library(tidyr)

# PLOT SETTINGS
tema=theme(axis.text.x = element_text(color="black",size=12, angle=0, family="Times New Roman",hjust=0.95,vjust=0.2),
           axis.text.y = element_text(color="black",size=12, family="Times New Roman"),
           axis.title = element_text(color="black",size=12, family="Times New Roman"),
           legend.text = element_text(color = "black",size=12, family="Times New Roman"),
           panel.border =element_rect(color = "black", fill = NA) , #element_blank(),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "left")
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

# OBTAIN ALL THE TAXA IN ALL THE RANKS FROM ALL THE METAGENOMES (1 TIME)
d=c("Bacteria","Archaea");p=c();c=c();o=c();f=c()

for (i in files){ print(keys[i]) #cycle by samples
  #load data
  full.taxa=fread(file = paste(data.dir,i,grep("phylodist",list.files(paste(data.dir,i,sep ="/" )),value = T),sep = "/"),sep = "\t")
  full.taxa=full.taxa[c(grep("Archaea",full.taxa$V5),grep("Bacteria",full.taxa$V5)),] #only prokaryotes
  full.taxa=cSplit(indt = full.taxa, splitCols = "V5", sep = ";", drop = TRUE)
  
  #get each taxa
  p=unique(c(p,unique(as.character(full.taxa$V5_2))))
  c=unique(c(c,unique(as.character(full.taxa$V5_3))))
  o=unique(c(o,unique(as.character(full.taxa$V5_4))))
  f=unique(c(f,unique(as.character(full.taxa$V5_5))))
  };remove(full.taxa)

save(d,file=paste(result.dir,"meta_domain_ids",sep="/"))
save(p,file=paste(result.dir,"meta_phylum_ids",sep="/"))
save(c,file=paste(result.dir,"meta_class_ids",sep="/"))
save(o,file=paste(result.dir,"meta_order_ids",sep="/"))
save(f,file=paste(result.dir,"meta_family_ids",sep="/")) #save all taxa

# COUNT EACH TAXA BASED ON THE CONTIG AVERAGE DEEPTH (MIGHT TAKE SOME TIME)
load(file=paste(result.dir,"meta_domain_ids",sep="/"));domain.data=data.frame(lineage=sort(d))
load(file=paste(result.dir,"meta_phylum_ids",sep="/"));phylum.data=data.frame(lineage=sort(p))
load(file=paste(result.dir,"meta_class_ids",sep="/"));class.data=data.frame(lineage=sort(c))
load(file=paste(result.dir,"meta_order_ids",sep="/"));order.data=data.frame(lineage=sort(o))
load(file=paste(result.dir,"meta_family_ids",sep="/"));family.data=data.frame(lineage=sort(f)) #load taxa

#thr=50 #select a threshold?

for (i in files){ print(keys[i])
  
  #load data
  full.taxa=fread(file = paste(data.dir,i,grep("phylodist",list.files(paste(data.dir,i,sep ="/" )),value = T),sep = "/"),sep = "\t",col.names = taxa.col)
  full.taxa=full.taxa[c(grep("Archaea",full.taxa$taxonomy),grep("Bacteria",full.taxa$taxonomy)),]
  full.depth=fread(file = paste(data.dir,i,grep("depth",list.files(paste(data.dir,i,sep ="/" )),value = T),sep = "/"), sep = "\t",col.names = depth.col)

  dep=full.depth$depth_value ; names(dep)=full.depth$orig.id #creates a vector with the depth and contig
  full.taxa$depth=dep[substr(full.taxa$gene_id, start = 1, stop = nchar(names(dep[1])))] #merges depth and annotation
  
  #if(is.numeric(thr)){full.taxa=full.taxa[full.taxa$percent_identity>thr,]} #similarity threhold

  d.count=c()
  for (j in sort(d)){ # Count each domain considering the depth of the contig
    print(j)
    d.count=c(d.count,sum(full.taxa[grep(j,full.taxa$taxonomy),depth])) #sums the dept of that gene in the contig
  } ; domain.data=cbind(domain.data,as.data.frame(d.count)) # Merges the result
  
  p.count=c()
  for (j in sort(p)){ # Count each phylum  considering the depth of the contig
    print(j)
    p.count=c(p.count,sum(full.taxa[grep(j,full.taxa$taxonomy),depth])) #sums the dept of that gene in the contig
  } ; phylum.data=cbind(phylum.data,as.data.frame(p.count)) # Merges the result
  
  c.count=c()
  for (j in sort(c)){ # Count each class  considering the depth of the contig
    print(j)
    c.count=c(c.count,sum(full.taxa[grep(j,full.taxa$taxonomy),depth])) #sums the dept of that gene in the contig
  } ; class.data=cbind(class.data,as.data.frame(c.count)) # Merges the result
  
  o.count=c()
  for (j in sort(o)){ # Count each order  considering the depth of the contig
    print(j)
    o.count=c(o.count,sum(full.taxa[grep(j,full.taxa$taxonomy),depth])) #sums the dept of that gene in the contig
  } ; order.data=cbind(order.data,as.data.frame(o.count)) # Merges the result
  
  f.count=c()
  for (j in sort(f)){ # Count each family  considering the depth of the contig
    print(j)
    f.count=c(f.count,sum(full.taxa[grep(j,full.taxa$taxonomy),depth])) #sums the dept of that gene in the contig
  } ; family.data=cbind(family.data,as.data.frame(f.count)) # Merges the result
  
}

colnames(domain.data)=c("lineage",keys[files])
colnames(phylum.data)=c("lineage",keys[files])
colnames(class.data)=c("lineage",keys[files])
colnames(order.data)=c("lineage",keys[files])
colnames(family.data)=c("lineage",keys[files]) #get column names

# SAVE DATA
write.table(domain.data, file = paste(result.dir,"domain_absolut_counts.txt",sep="/"),sep="\t", col.names = T, row.names = F)
write.table(phylum.data, file = paste(result.dir,"phylum_absolut_counts.txt",sep="/"),sep="\t", col.names = T, row.names = F)
write.table(class.data, file = paste(result.dir,"class_absolut_counts.txt",sep="/"),sep="\t", col.names = T, row.names = F)
write.table(order.data, file = paste(result.dir,"order_absolut_counts.txt",sep="/"),sep="\t", col.names = T, row.names = F)
write.table(family.data, file = paste(result.dir,"family_absolut_counts.txt",sep="/"),sep="\t", col.names = T, row.names = F)

