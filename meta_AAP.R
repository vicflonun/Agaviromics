#CREATES ABUNDANCE PLOTS FOR AEROBIC ANOXYGENIC PHOTOTROPHY GENES #

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

#all
guild=list(name="AAP",genes=c("bch","chl","puf","crt","E1.14.13.81","puc","puh"))
function_id=data.frame()
for (i in guild$genes) {function_id=unique(rbind(function_id,path.data[grep(i,path.data$anygene_name),c("ko_id","anygene_name")] ))}
ids=function_id$anygene_name ; names(ids)=function_id$ko_id

#by process
markers=list(Reaction_center=c("pufM", "pufL", "pufA" ,"pufB","puhA"),
             Chlorophyll=c("chlG","chlP"),
             Carotenes=c("crtD","crtF","crtC"))

function_id=data.frame()
for (i in names(markers)) {function_id=rbind(function_id,path.data[path.data$anygene_name %in% markers[[i]],c("ko_id","anygene_name")] )}
function_id=unique(function_id)


# CRATES A HEATMAP OF GENE ABUNDANCES
library(pheatmap); library(colorspace)

# extract gene markers
result_norm=as.data.frame(t(t(data.full)/colSums(data.full)))
mark=data.frame(row.names = colnames(data.full))
for (j in names(markers)){
  nombres = ids[ids %in% markers[[j]]]
  sub=result_norm[rownames(result_norm) %in% names(nombres),]
  rownames(sub)=nombres[rownames(sub)]
  mark=cbind(mark,t(sub*1e6))
  
}

# ORDER ROWS AND COLS
heat=t(mark[rev(orden),])
orden2=c(markers$Reaction_center,markers$Chlorophyll,markers$Carotenes)
orden2=orden2[orden2 %in% rownames(heat)]

# PLOT
png(filename = paste(result.dir,"AAP_HEAT_2.png",sep = "/"), 
    width = 100, height = 100, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")

pheatmap(t(heat[orden2,]), color = c("white",sequential_hcl(palette = "grays",n=3000,rev = T)), 
         cluster_rows = F, cluster_cols = F,
         cellwidth = 13, cellheight = 9, border_color = "white")
dev.off()



# COUNT TAXA OF SELECTED GENES

# count
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
write.table(result,paste(result.dir,paste(guild$name,"taxa_by_function.txt",sep = "_"), sep = "/"), sep = "\t", row.names = F, col.names = T, quote = F)

# PLOT TAXA COUNTS
library(tidyr);library(Rmisc)

# load data
result=read.delim(file=paste(result.dir,paste(guild$name,"taxa_by_function.txt",sep = "_"), sep = "/"))
result=result[result$func_id %in% function_id$ko_id,]
lineages=sort(as.character(unique(result$order)))

# transform data for plotting
data=data.frame(row.names = lineages)
for (i in keys[files]){
  for (j in lineages){data[j,i]=sum(result[result$sample == i & result$order==j, "depth"])}
}

result_norm=as.data.frame(t(t(data)/colSums(data)*100)) #relative abundance
result_norm$lineage=rownames(result_norm)
result.tidy=gather(data=result_norm, key = assembled.library, value = relative.abundance, keys[files])
result.tidy[is.nan(result.tidy$relative.abundance),"relative.abundance"]=0

# order levels
lin.orden=summarySE(result.tidy, measurevar = "relative.abundance", groupvars = c("lineage")) #mean abundance in all metagenomes
lin.orden=as.character(lin.orden[order(lin.orden$relative.abundance,decreasing = T),"lineage"]) #sort

#filter by abundance
filter=5
result.tidy$lineage=factor(result.tidy$lineage, levels = c(lin.orden,"Other",paste("Lineages<",filter,"%",sep = ""))) #set levels with the filter
result.tidy[result.tidy$relative.abundance<filter,"lineage"]=paste("Lineages<",filter,"%",sep = "") #filter
imp=grep("[:.:]e",unique(result.tidy$assembled.library), value = T)
taximp=as.character(unique(result.tidy[result.tidy$assembled.library %in% imp, "lineage"]))
result.tidy[!(result.tidy$lineage %in% taximp),"lineage"]="Other"

#color
ncolor = length(unique(result.tidy$lineage));print(ncolor) #sanity check for the numer of taxa (better less than 12)

#plot
png(filename = paste(result.dir,"AAP_barras_3.png",sep = "/"), 
    width = 100, height = 80, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")

ggplot(data=result.tidy, aes(x=assembled.library, y=relative.abundance, fill=lineage))+
  geom_bar(stat = "identity")+scale_fill_manual(values = c(my_pal(ncolor)[1:(ncolor-2)],"grey60","black"))+
  scale_x_discrete(limits=orden)+tema+ylab("Relative abundance (%)")+coord_flip()

dev.off()







