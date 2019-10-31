### COUNT THE GENES AND PATHWAYS WITH KO ANNOTATION IN EACH METAGENOME BY LINEAGE ##

# WORKING DIRECTORIES
setwd("C:/my_path")
data.dir=c("C:/my_data") #metagenomic data
result.dir=c("C:/my_path/counts_per_taxa_domain")
dir.create(result.dir, recursive = T, showWarnings = F)

# LIBRARIES
library(data.table)

# FUNCTION 
func="ko"

# COLUMN NAMES OF DATA
taxa.col=c("gene_id","homolog_gene_oid", "homolog_taxon_oid", "percent_identity","taxonomy")
if (func=="ko"){anot.col=c("gene_id","img_ko_flag", "func_id" ,"percent_identity" ,"query_start","query_end","subj_start","subj_end","evalue","bit_score","align_length")}
depth.col=c("orig.id","depth")

# LOAD METADATA
path.data=fread(file = paste(getwd(),"kegg_pathway_data.txt",sep="/"),sep="\t", header = T) #pathway data
paths=sort(unique(path.data$pathway)) #all pathways
process.data=fread(file = paste(getwd(),"kegg_process_data.txt",sep="/"),sep="\t", header = T) #process data
full.meta=read.delim(paste(data.dir,"meta_data.txt",sep = "/"),header = T,stringsAsFactors = F) #metagenome data
load(file=paste(getwd(),"ko.ids",sep="/")) #ko data
ids=gsub("KO:","",ids) #modify ko

# GET FILE NAMES 
files=grep(33000,list.files(data.dir),value = T)
keys=full.meta$Genome_name ; names(keys)=full.meta$taxon_oid

 
# LOAD THE NECESSARY TAXA
load(file=paste(getwd(),"meta_domain_ids",sep="/"))
load(file=paste(getwd(),"meta_phylum_ids",sep="/"))

# SELECT THE TAXA FOR COUNTING (USE YOUR OWN COMMANDS)
taxon=c(p)

# COUNT #
# result will count each ko category of each taxa in each metagenome (list of 20 tables)
# p.result will count each pathway of each taxa in each metagenome (list of 20 tables)
result=list()
p.result=list()

for (i in files){ print(keys[i])
  
  #read data
  full.anot=fread(file = paste(data.dir,i,grep(func,list.files(paste(data.dir,i,sep ="/" )),value = T),sep = "/"),sep = "\t", col.names = anot.col)
  full.anot$func_id=gsub("KO:","",full.anot$func_id)  #change the names of the func_id
  full.depth=fread(file = paste(data.dir,i,grep("depth",list.files(paste(data.dir,i,sep ="/" )),value = T),sep = "/"), sep = "\t",col.names = depth.col)
  full.taxa=fread(file = paste(data.dir,i,grep("phylodist",list.files(paste(data.dir,i,sep ="/" )),value = T),sep = "/"),sep = "\t",col.names = taxa.col)
  dep=c(full.depth$depth) ; names(dep)=full.depth$orig.id #create a vector with the depth and contig
  full.anot$depth=dep[substr(full.anot$gene_id, start = 1, stop = nchar(names(dep[1])))] #merge depth and annotation
  
  #generate counting tables
  count.tab=data.frame(row.names = sort(ids))
  p.count.tab=data.frame(row.names = sort(paths))
  
  for (k in taxon){ print(k) #for each taxa
    
    #extract functions of specific taxa
    extract=full.taxa[grep(k,full.taxa$taxonomy),gene_id]
    sub.anot=full.anot[full.anot$gene_id %in% extract,]
    
    count=c()
    for (j in sort(ids)){
      count[j]=sum(sub.anot[sub.anot$func_id==j,depth])
      } #count each KO category 
    
    p.count=c()
    for (j in sort(paths)){
      extract2=path.data[path.data$pathway==j,ko_id]
      p.count[j]=sum(sub.anot[sub.anot$func_id %in% extract2,depth])
      }  #count each paths category 
    
    count.tab=data.frame(count.tab,cuenta=count)  
    p.count.tab=data.frame(p.count.tab,cuenta=p.count) #merge the counts of each taxa
  }
    
    #name each colum with its especific taxa
    colnames(count.tab)=taxon 
    colnames(p.count.tab)=taxon
    
    #merges the table of each metagenome
    result[[i]]=count.tab 
    p.result[[i]]=p.count.tab
}

#save
save(result, file=paste(result.dir,"indirect.taxon.ko.counts.list",sep = "/"))
save(p.result, file=paste(result.dir,"indirect.taxon.path.counts.list",sep = "/"))

# TRANSFORM
#Tranform the data sets from tables of metagenomes to tables of taxas (n taxa tables with 20 columns)

#KO categories
taxon.tab.list=list()
for (j in 1:length(taxon)){
  taxon.tab=data.frame(row.names = sort(ids))
  for ( i in files){
    taxon.tab=data.frame(taxon.tab,count=result[[i]][,j])
  }
  colnames(taxon.tab)=keys[files]
  taxon.tab.list[[taxon[j]]]=taxon.tab
}
save(taxon.tab.list, file=paste(result.dir,"taxon.ko.counts.list",sep = "/"))

#pathways
p.taxon.tab.list=list()
for (j in 1:length(taxon)){
  p.taxon.tab=data.frame(row.names = sort(paths))
  for ( i in files){
    p.taxon.tab=data.frame(p.taxon.tab,count=p.result[[i]][,j])
  }
  colnames(p.taxon.tab)=keys[files]
  p.taxon.tab.list[[taxon[j]]]=p.taxon.tab
}

save(p.taxon.tab.list, file=paste(result.dir,"taxon.paths.counts.list",sep = "/"))

