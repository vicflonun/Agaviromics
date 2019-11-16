# CREATES ABUNDANCE PLOTS FOR BIOFILM FORMATION AND QS GENES #

# WORKING DIRECTORIES
setwd("/my_dir/")
data.dir=c("my_data") #metagenomic data
result.dir=c("/my_dir//result")
dir.create(result.dir, recursive = T, showWarnings = F)

# LIBRARIES
library(data.table)
library(MASS)
library(vegan)
library(ggplot2)
library(ggrepel)
library(extrafont)
#font_import()
loadfonts(device = "win")
names(wf[wf=="TT Times New Roman"])


# plot settings
tema=theme(axis.text.x = element_text(color="black",size=12, angle=0, family="Times New Roman"),#,hjust=0.95,vjust=0.2),
           axis.text.y = element_text(color="black",size=12, family="Times New Roman"),
           axis.title = element_text(color="black",size=12, family="Times New Roman"),
           legend.text = element_text(color = "black",size=12, family="Times New Roman"),
           panel.border =element_rect(color = "white", fill = NA) , #element_blank(),
           panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right")

# set colors 
colo=c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","grey20")
my_pal = colorRampPalette(colo)

# set sample order
orden=c("Am.s","Pe.s","Ma.s","Sf.s","At.Am.rz","At.Pe.rz","As.Ma.rz","As.Sf.rz", "Mg.Ma.rz","Mg.Sf.rz" ,
        "Or.Ma.rz" ,"Or.Sf.rz",  "At.Am.e","At.Pe.e","As.Ma.e","As.Sf.e","Mg.Ma.e","Mg.Sf.e", "Or.Ma.e" ,"Or.Sf.e") 


#ARGUMENTS
func="ko"

#SET LEVELS AND COLUMN NAMES OF DATA
taxa.col=c("gene_id","homolog_gene_oid", "homolog_taxon_oid", "percent_identity","taxonomy")
if (func=="ko"){anot.col=c("gene_id","img_ko_flag", "func_id" ,"percent_identity" ,"query_start","query_end","subj_start","subj_end","evalue","bit_score","align_length")}
depth.col=c("orig.id","depth")

path.data=read.delim(file = paste(getwd(),"kegg_pathway_data.txt",sep="/"),sep="\t", header = T, quote = "",stringsAsFactors = F)
full.meta=read.delim(paste(data.dir,"meta_data.txt",sep = "/"),header = T,stringsAsFactors = F)

# GET FILE NAMES 
files=grep(33000,list.files(data.dir),value = T)
keys=full.meta$Genome_name ; names(keys)=full.meta$taxon_oid

# LOAD COUNT DATA
load(file=paste(getwd(),"counts_per_taxa_domain","taxon.ko.counts.list",sep = "/"))
data.full=taxon.tab.list[["Bacteria"]]+taxon.tab.list[["Archaea"]]

# SUBSET GENES

# all
bio_genes=path.data[path.data$pathway %in% c("Quorum sensing [PATH:ko02024]","Biofilm formation - Pseudomonas aeruginosa [PATH:ko02025]",
          "Biofilm formation - Escherichia coli [PATH:ko02026]","Flagellar assembly [PATH:ko02040]",
          "Biofilm formation - Vibrio cholerae [PATH:ko05111]",
          "Bacterial chemotaxis [PATH:ko02030]"),c("ko_id", "anygene_name")]

guild=list(name="BIO.QS",genes= bio_genes$anygene_name)

# or only biofilm enriched genes from meta_enrichment.R
phyllo=read.delim(file="C:/Users/victo/Google Drive/cam_final_comp_2019/PAPER/SYM_Prokaryotes_ko_75/phyllo_vs_soil/Expanded_families_edgeR.txt")
rhizo=read.delim(file="C:/Users/victo/Google Drive/cam_final_comp_2019/PAPER/SYM_Prokaryotes_ko_75/rhizo_vs_soil/Expanded_families_edgeR.txt")
p=phyllo[rownames(phyllo) %in% bio_genes$ko_id,]
r=rhizo[rownames(rhizo) %in% bio_genes$ko_id,]

guild=list(name="BIO.QS",genes=unique(path.data[path.data$ko_id %in% c(rownames(p), rownames(r)),"anygene_name"]))

function_id=data.frame()
for (i in guild$genes) {function_id=unique(rbind(function_id,path.data[grep(i,path.data$anygene_name),c("ko_id","anygene_name")] ))}
ids=function_id$anygene_name ; names(ids)=function_id$ko_id

## CRATES A HEATMAP OF GENE ABUNDANCES
library(pheatmap); library(colorspace)

# select genes
markers=c(c("luxS","lsrR","lasI"),
          c("pgaA","rfbF","srfAB"),
          c("arcA","flhC","rcsB"))
# get ids of genes
function_id=data.frame()
for (i in markers) {function_id=unique(rbind(function_id,path.data[grep(i,path.data$anygene_name),c("ko_id","anygene_name")] ))}
ids=function_id$anygene_name ; names(ids)=function_id$ko_id
ids=ids[!names(ids) %in% c("K01478","K11029","K22944","K00978") ]

# subset counts
result_norm=as.data.frame(t(t(data.full)/colSums(data.full)))
mark=data.frame(row.names = colnames(data.full))
for (j in markers){
  nombres = ids[ids %in% j]
  sub=result_norm[rownames(result_norm) %in% names(nombres),]
  rownames(sub)=nombres[rownames(sub)]
  mark=cbind(mark,t(sub*1e6))
}

# order cols and rows
heat=t(mark[rev(orden),])

library(pheatmap);library(colorspace)
#save plot
png(filename = paste(result.dir,"BIO_HEAT_3.png",sep = "/"), 
    width = 150, height = 100, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")

pheatmap(t(heat), color = c("white",sequential_hcl(palette = "grays",n=3000,rev = T)), 
         cluster_rows = F, cluster_cols = F,
         cellwidth = 13, cellheight = 9, border_color = "white")

dev.off()


# COUNT TAXA OF SELECTED GENES

result=data.frame()
for (i in files){ print(keys[i])
  
  # read data
  full.anot=fread(file = paste(data.dir,i,grep(func,list.files(paste(data.dir,i,sep ="/" )),value = T),sep = "/"),sep = "\t", col.names = anot.col)[,c(1,3)]
  full.anot$func_id=gsub("KO:","",full.anot$func_id)  #change the names of the func_id
  full.depth=fread(file = paste(data.dir,i,grep("depth",list.files(paste(data.dir,i,sep ="/" )),value = T),sep = "/"), sep = "\t",col.names = depth.col)
  full.taxa=fread(file = paste(data.dir,i,grep("phylodist",list.files(paste(data.dir,i,sep ="/" )),value = T),sep = "/"),sep = "\t",col.names = taxa.col)[,c(1,5)]
  dep=c(full.depth$depth) ; names(dep)=full.depth$orig.id #Creates a vector with the depth and contig
  full.anot$depth=dep[substr(full.anot$gene_id, start = 1, stop = nchar(names(dep[1])))] #merges depth and annotatio
  
  # subset data
  full.taxa=full.taxa[c(grep("Archaea",full.taxa$taxonomy),grep("Bacteria",full.taxa$taxonomy)),]
  full.anot=full.anot[full.anot$gene_id %in% full.taxa$gene_id,]
  anot=full.anot[full.anot$func_id %in% function_id$ko_id,]
  taxa=full.taxa[full.taxa$gene_id %in% anot$gene_id,]
  
  if(dim(taxa)[1]>0){
  # replace data
  sub.result=data.frame(anot, taxonomy=rep(paste("NT","NT","NT","NT","NT","NT","NT","NT",sep = ";")), stringsAsFactors = F)
  for (j in taxa$gene_id){sub.result[grep(j,sub.result$gene_id),"taxonomy"]=taxa[taxa$gene_id==j,taxonomy]}
  # add metadata
  meta=data.frame(sample=rep(keys[i], times=dim(sub.result)[1]),
                  species=rep(full.meta[full.meta$taxon_oid==i,"Species"], times=dim(sub.result)[1]),
                  compartment=rep(full.meta[full.meta$taxon_oid==i,"Compartment"], times=dim(sub.result)[1]),
                  site=rep(full.meta[full.meta$taxon_oid==i,"Site"], times=dim(sub.result)[1]))
  #merge
  result = rbind(result, cbind(sub.result,meta)) } else {print(paste(keys[i],"function not found"))}
  
  }
library(splitstackshape)
result=cSplit(indt = result, splitCols = "taxonomy", sep = ";", drop = TRUE)
colnames(result)[8:15]=c("domain","phylum","class","order","family","genus","species","strain")
write.table(result,paste(result.dir,paste(guild$name,"marker_taxa_by_function.txt",sep = "_"), sep = "/"), sep = "\t", row.names = F, col.names = T, quote = F)


# PLOT TAXA COUNTS
library(tidyr)

#markers
markers=c(c("luxS","lsrR","lasI"),
          c("pgaA","rfbF","srfAB"),
          c("arcA","flhC","rcsB"))
#load data
result=read.delim(file=paste(result.dir,paste("biofilm","marker_taxa_by_function.txt",sep = "_"), sep = "/"))

#subset data
mark=unlist(markers)
result=result[result$func_id %in% path.data[path.data$anygene_name %in% mark, "ko_id"],]
function_id=data.frame()
for (i in mark) {function_id=rbind(function_id,path.data[path.data$anygene_name==i,c("ko_id","anygene_name")] )}
function_id=unique(function_id) ; print(function_id)
function_id=function_id[!function_id$ko_id %in% c("K01478","K01478","K11029","K22944","K00978","K00978"),]

# organize data for plotting
lineages=sort(as.character(unique(result$order)))
result.full=data.frame()
for (l in function_id$ko_id){
  sub.result=result[result$func_id == l,]  
  sub.data=data.frame(row.names = lineages)
  for (i in keys[files]){
    for (j in lineages){sub.data[j,i]=sum(sub.result[sub.result$sample == i & sub.result$order==j, "depth"])}
  }
  result_norm=as.data.frame(t(t(sub.data)/colSums(sub.data)))
  #result_norm=sub.data
  result_norm$lineage=rownames(result_norm)
  result_norm$marker=rep(function_id[function_id$ko_id==l,"anygene_name"])
  result_norm=gather(data=result_norm, key = "sample", value = "relative.abundance", keys[files])
  result.full=rbind(result.full,result_norm)
}

result.tidy=result.full[result.full$marker %in% mark,]
#result.tidy=result.full
result.tidy[is.nan(result.tidy$relative.abundance),"relative.abundance"]=0

#order levels
lin.orden=summarySE(result.tidy, measurevar = "relative.abundance", groupvars = c("lineage")) #mean abundance in all metagenomes
lin.orden=as.character(lin.orden[order(lin.orden$relative.abundance,decreasing = T),"lineage"]) #sort

#filter by abundance
filter=0.3
result.tidy$lineage=factor(result.tidy$lineage, levels = c(lin.orden,"Other",paste("Lineages<",filter,"%",sep = ""))) #set levels with the filter
result.tidy[result.tidy$relative.abundance<filter,"lineage"]=paste("Lineages<",filter,"%",sep = "") #filter
imp=c(grep("[:.:]e",unique(result.tidy$sample), value = T),grep("[:.:]rz",unique(result.tidy$sample), value = T))
taximp=as.character(unique(result.tidy[result.tidy$sample %in% imp, "lineage"]))
result.tidy[!(result.tidy$lineage %in% taximp),"lineage"]="Other"

#color
ncolor = length(unique(result.tidy$lineage));print(ncolor) #sanity check for the numer of taxa (better less than 12)

result.tidy$marker=factor(result.tidy$marker, levels = markers)


#save plot
png(filename = paste(result.dir,"BIO_barras_6.png",sep = "/"), 
    width = 200, height = 90, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")

ggplot(data=result.tidy, aes(x=sample, y=relative.abundance*100, fill=lineage))+
  geom_bar(stat = "identity")+scale_fill_manual(values = my_pal(ncolor))+
  scale_x_discrete(limits=orden)+
  scale_y_continuous(breaks = c(25,75))+
  tema+ylab("Relative abundance (%)")+
  facet_wrap(~marker)+
  coord_flip()

dev.off()

