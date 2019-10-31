# ENRICHMENT ANALYSIS OF KO CATEGORIES #

# LIBRARY
library(edgeR)

# SETTINGS

setwd(paste("C:/my_path"))
full.meta=read.delim(paste("C:/my_data","meta_data.txt",sep = "/"),header = T,stringsAsFactors = F)

logFC_thres = 0    #folding change threshold
FDR_thres = 0.05   #FDR threshold
funcion="ko"        #function
type="RHI"          #comparison
                    #ALL = comparrison btw compartments of all plants
                    #SYM = comparrison btw compartments of all sympatric plants
                    #TEQ = comparrison btw compartments of A. tequilana
                    #SAL = comparrison btw compartments of A. salmiana
                    #PHY = comparrison btw species of phyllospheres
                    #RHI = comparrison btw species of rhizospheres
cero=c(0.25,0.50,0.75) #prevalence threshold in percentage, can use multiple
folder="TEST"


#taxa selection
taxa=c("Prokaryotes","Cyanobacteria","Firmicutes", "Proteobacteria","Bacteroidetes","Actinobacteria")[1]

for (samcero in cero){ #for each prevalence threshols
 
  # DATA   
   
#set directories
main_dir=paste(getwd(),folder,paste(type,taxa,funcion,samcero,sep = "_"),sep = "/") #to save result
dir.create(main_dir, recursive = T, showWarnings = F)
desc=paste(funcion,"ids","edgeR.txt",sep= "_") #ko description names

#load data
if (taxa=="Prokaryotes"){
  data.dir=c("C:/gene_couns_path/counts_per_taxa_domain")
  load(file=paste(data.dir,"taxon.ko.counts.list",sep = "/"))
  count_tab_full=taxon.tab.list[["Bacteria"]]+taxon.tab.list[["Archaea"]]
} else {
  data.dir = c("C:/gene_couns_path//counts_per_taxa_phylum")
  load(file=paste(data.dir,"taxon.ko.counts.list",sep = "/"))
  count_tab_full=taxon.tab.list[[taxa]]
}

  # ORGANIZE SAMPLES 
if (type=="SYM"){
  group_1 = c("As.Ma.e","As.Sf.e","Mg.Ma.e","Mg.Sf.e","Or.Sf.e","Or.Ma.e")
  group_2 = c("As.Ma.rz","As.Sf.rz","Mg.Ma.rz","Mg.Sf.rz","Or.Ma.rz","Or.Sf.rz")
  group_3 = c("Ma.s","Sf.s" )
} else if (type=="TEQ"){
  group_1 = c("At.Am.e","At.Pe.e") ; group_2 = c("At.Am.rz","At.Pe.rz") ; group_3 = c("Am.s","Pe.s" )
} else if (type=="PHY"){
  group_1 = c("As.Ma.e","As.Sf.e") ; group_2 = c("Mg.Ma.e","Mg.Sf.e","Or.Sf.e","Or.Ma.e")
  group_3 = c("At.Am.e","At.Pe.e")
} else if (type=="ALL"){
  group_1 = c("As.Ma.e","As.Sf.e","Mg.Ma.e","Mg.Sf.e","Or.Sf.e","Or.Ma.e","At.Am.e","At.Pe.e")
  group_2 = c("As.Ma.rz","As.Sf.rz","Mg.Ma.rz","Mg.Sf.rz","Or.Ma.rz","Or.Sf.rz","At.Am.rz","At.Pe.rz")
  group_3 = c("Ma.s","Sf.s","Am.s","Pe.s")
} else if (type=="SAL"){
  group_1 = c("As.Ma.e","As.Sf.e") ; group_2 = c("As.Ma.rz","As.Sf.rz") ; group_3 = c("Ma.s","Sf.s")
} else if (type=="RHI"){
  group_1 = c("As.Ma.rz","As.Sf.rz") ; group_2 = c("Mg.Ma.rz","Mg.Sf.rz","Or.Sf.rz","Or.Ma.rz")
  group_3 = c("At.Am.rz","At.Pe.rz")
}

  # ORDER GROUPS
agave_order = c(group_1, group_2,group_3)
count_tab_full = count_tab_full[,agave_order]

  # FILTER BY ABUNDANCE
filter=mean(rep(1,times=length(agave_order))/colSums(count_tab_full)*1e6) # calculates how many cpm is a read
count_tab = count_tab_full[rowSums(cpm(count_tab_full) > filter) >= length(agave_order)*samcero,] #filter
#count_tab = count_tab_full #no filter

  # LOAD DESCRIPTION 
annot_tab=read.delim(paste(getwd(),desc,sep = "/"))
rownames(annot_tab)=annot_tab$func_id
colnames(annot_tab)=c("func_id" ,"desc" , "desc_id")
  
  # REMOVE UNCHARACTERIZED PROTEINS
if (funcion=="cog"){
  removing=as.character(annot_tab[grep("Uncharacterized",annot_tab$desc),"func_id"])
  count_tab=count_tab[!(rownames(count_tab)%in%removing),]
} else { 
  removing=as.character(annot_tab[grep("Other",annot_tab$desc_id),"func_id"])
  count_tab=count_tab[!(rownames(count_tab)%in%removing),]
}

  # SET THE COMPARISONS
if (type=="TEQ" | type=="SYM" | type=="ALL"|type=="SAL"){
  factors_setting = c(rep("phyll",times = length(group_1)), 
                     rep("rhizo",times = length(group_2)),
                     rep("soil",times = length(group_3))) #set a level for each sample
  grp = factor(factors_setting)  #set as factors 
  design = model.matrix(~0+grp) #create a desing matrix
  colnames(design) = levels(grp) #change colnames
  contrasts = makeContrasts( 
    "phyllo_vs_rhizo" = "(phyll)/1 - (rhizo)/1",
    "phyllo_vs_soil" = "(phyll)/1 - (soil)/1",
    "rhizo_vs_soil" = "(rhizo)/1 - (soil)/1",levels=design) #define contrasts of interest

} else if (type=="PHY" | type=="RHI"){
  factors_setting = c(rep("salmiana",times = length(group_1)), 
                       rep("cacti",times = length(group_2)),
                       rep("tequilana",times = length(group_3))) #set a level for each sample
  grp = factor(factors_setting) #set as factors 
  design = model.matrix(~0+grp) #create a desing matrix
  colnames(design) = levels(grp) #change colnames
  contrasts = makeContrasts(
    "salmiana_vs_tequilana" = "(salmiana)/1 - (tequilana)/1",
    "salmiana_vs_cacti" = "(salmiana)/1 - (cacti)/1",
    "cacti_vs_tequilana" = "(cacti)/1 - (tequilana)/1",levels=design) #define contrasts of interest
}

# CREATE THE EDGER OBJECT TO WORK WITH
cg = DGEList(counts=count_tab, group=grp)
cg = calcNormFactors(cg) #normalize counts
#plotMDS.DGEList(cg)
cg = estimateGLMCommonDisp(cg, design, verbose=T) #estimate dispersion 
venn_tab = data.frame(row.names = rownames(cg$counts)) #result table
fit = glmFit(cg, design) #fit to the statistical model

# COMPARISONS 

for (contrast in colnames(contrasts)) {
  
  print(contrast)
  contrast_dir = paste(main_dir, contrast, sep = "/")
  dir.create(contrast_dir, showWarnings = F) #create output directory
  
  lrt = glmLRT(fit, contrast=contrasts[,contrast]) #performs the likelyhood ratio tests
  venn_col = as.data.frame(decideTestsDGE(lrt, adjust.method = "BH", p.value = FDR_thres)) #indentify the differential genes in a matrix
  colnames(venn_col) = contrast
  venn_tab = cbind(venn_tab, venn_col) #merge dataframes
  cgTab = topTags(lrt, n=Inf)$table #select top differential genes
  cgTab = signif(cgTab, digits = 3) #reduce decimals
  cgTab$desc = annot_tab[rownames(cgTab),]$desc #add description 
  cgTab$desc_id = annot_tab[rownames(cgTab),]$desc_id #add description 
  deGenes = rownames(cgTab)[cgTab$FDR < FDR_thres & abs(cgTab$logFC) > logFC_thres] #select the differential genes based on threshold
  poscgTab = cgTab[(cgTab$FDR < FDR_thres & abs(cgTab$logFC) >= logFC_thres & cgTab$logFC>0),] #enriched
  negcgTab = cgTab[(cgTab$FDR < FDR_thres & abs(cgTab$logFC) >= logFC_thres & cgTab$logFC<0),] #reduced
  
  #Print result
  posDE = nrow(poscgTab)
  print("expanded:")
  print(posDE)
  negDE = nrow(negcgTab)
  print("reduced:")
  print(negDE)
  
  #merge counts and statistical information
  poscgTab = cbind(poscgTab, cg$counts[rownames(poscgTab),])
  negcgTab = cbind(negcgTab, cg$counts[rownames(negcgTab),])
  
  #create threshold directories
  thres_dir = contrast_dir
  dir.create(thres_dir, showWarnings = F)
  
  #save
  write.table(cgTab, file=paste(thres_dir, "Comparative_Genomics_families_edgeR.txt", sep = "/"), sep="\t", quote=F)
  write.table(poscgTab, file=paste(thres_dir, "Expanded_families_edgeR.txt", sep = "/"), sep="\t", quote=F)
  write.table(negcgTab, file=paste(thres_dir, "Reduced_families_edgeR.txt", sep = "/"), sep="\t", quote=F)
}

}


# CREATES EASY TO READ FILES FROM ENRICHMENT ANALYSIS

# LIBRARIES
library(ggplot2)
library(splitstackshape)

# DIRECTORIES
setwd("C:/my_path")
main.dir="C:/my_path/result"
path.data=read.delim(paste(getwd(),"new_kegg_pathway_data.txt",sep = "/"), quote = "")

# FILES TO USE
comp=grep("ko",list.files(main.dir),value = T)

# FILES NAMES
files=c("Expanded_families_edgeR.txt","Reduced_families_edgeR.txt")

final.result=path.data
for (comparison in comp){ 
  print(comparison)
  
  if (sum(grep("phyllo",list.files(paste(main.dir,comparison,sep = "/"))))>0){
    pair=c("phyllo_vs_rhizo","phyllo_vs_soil", "rhizo_vs_soil")
  } else {
    pair=c("salmiana_vs_tequilana","salmiana_vs_cacti", "cacti_vs_tequilana")
  }
  
  result=path.data
  for (pairs in pair){
    print(pairs)
    name=strsplit(pairs,"_vs_",fixed = T)[[1]]
    
    sub.result=path.data
    sub.result[,5]="0"
    for (i in 1:2){
      
      print(name[i])
      data=read.delim(paste(main.dir,comparison,pairs,paste(files[i],sep = "."),sep = "/"), stringsAsFactors = F)
      sub.result[sub.result$ko_id %in% rownames(data),5]=name[i]
      
    }
    
    result=cbind(result,sub.result[,5])
    
  }
  colnames(result)[c(5:7)]=paste(comparison,pair,sep = "_")
  result=result[,c(5:7)]
  final.result=data.frame(final.result,result)
  print("________________________________________________________________________")
}

write.table(final.result, file = paste(main.dir,"meta.sumary.txt", sep = "/"), col.names = T, row.names = F, quote = F, sep = "\t")




