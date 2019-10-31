## OVERREPRESENTATION ANALYSIS ##

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


path.data=read.delim(paste(getwd(),"new_kegg_pathway_data.txt",sep = "/"), quote = "")
paths=unique(as.character(path.data$pathway))

#select comparison
for (comparison in comp){ print(comparison)
  
  if (sum(grep("phyllo",list.files(paste(main.dir,comparison,sep = "/"))))>0){
    pair=c("phyllo_vs_rhizo","phyllo_vs_soil", "rhizo_vs_soil")}
  else if ((sum(grep("native",list.files(paste(main.dir,comparison,sep = "/"))))>0)){
    pair=c("native_vs_cultivated")
  } else {
    pair=c("salmiana_vs_tequilana","salmiana_vs_cacti", "cacti_vs_tequilana")
  }
  
  #select each element of the comparison
  for (pairs in pair){
    print(pairs)
    name=strsplit(pairs,"_vs_",fixed = T)[[1]]
    
    if (file.exists(paste(main.dir,comparison,pairs,sep = "/"))){
      
      pvalor=data.frame(row.names = paths) #result dataframe
      for (i in 1:length(files)){
        
        #read data
        top=read.delim(paste(main.dir,comparison,pairs,files[i],sep = "/")) #differentially enriched
        all=read.delim(paste(main.dir,comparison,pairs,"Comparative_Genomics_families_edgeR.txt",sep = "/")) #all data
        
        pathways.top=data.frame(row.names=paths) #for enriched genes (selection)
        pathways.all=data.frame(row.names=paths) #for all genes (urn)
        pathways.prop=data.frame(row.names=paths) #for proportion data
        
        #count tne number of enzymes and in each pathway in the enriched genes
        for (j in 1:length(paths)){
          count=sum(as.character(top$desc_id) %in% unique(as.character(path.data[path.data$pathway==paths[j],"enzyme"])))
          pathways.top[j,"count"]=count
        }
        
        #count tne number of enzymes and its proportion in each pathway in all the genes 
        for (j in 1:length(paths)){
          count=sum(as.character(all$desc_id) %in% unique(as.character(path.data[path.data$pathway==paths[j],"enzyme"])))
          pathways.all[j,"count"]=count
          ids=as.character(path.data[path.data$pathway==paths[j],"ko_id"])
          prop=length(rownames(all[rownames(all) %in% ids ,]))/length(ids)*100
          pathways.prop[j,"proportion"]=prop
        }
        
        #performs the hypergeometric test
        p=c()
        for (j in paths) {
          p[j]=phyper(q=pathways.top[j,"count"] , #white balls drawn
                      m=pathways.all[j, "count"] , #white balls in the urn
                      n=dim(all)[1]-pathways.all[j, "count"] , #black ball in the urn
                      k=dim(top)[1] , #number of balls drawn
                      lower.tail = F, log.p = F)
        }
        p=p.adjust(p,method="BH",n=length(p)) #multiple testing correction
        p[which(pathways.top<=2)]=1 #pathways with very few genes not significant
        p[which(pathways.prop<=25)]=1 #pathways with low proportion of enzymes not significant
        p[which(p==0)]=1e-300 #correct that R can not oprate les than -300
        pvalor=data.frame(pvalor,p)
        
        print(name[i])
        print(as.data.frame(p[which(p<=0.05)])) #print result
      }
      
      #save
      colnames(pvalor)=name
      write.table(cbind(pathways.all,pathways.prop), file = paste(main.dir,comparison,pairs,paste("enriched.enzymes","txt",sep = "."),sep = "/"), quote = F, col.names = T, row.names = T, sep = "\t")
      write.table(pvalor, file = paste(main.dir,comparison,pairs,paste("OVA","txt",sep = "."),sep = "/"), quote = F, col.names = T, row.names = T, sep = "\t")
    } else {print("file is cero")}
  }
}

## SUMMARY OF THE OVA ANLYSIS ##

# FILES TO USE
comp=grep("ko",list.files(main.dir),value = T)
comp=grep("75",comp, value = T)[-c(1,3,4)]

ova.result=data.frame(row.names = sort(unique(path.data$pathway)))
for (comparison in comp){ print(comparison)
  
  #select comparisons
  if (sum(grep("phyllo",list.files(paste(main.dir,comparison,sep = "/"))))>0){
    pair=c("phyllo_vs_rhizo","phyllo_vs_soil", "rhizo_vs_soil")
    } else if ((sum(grep("native",list.files(paste(main.dir,comparison,sep = "/"))))>0)){
    pair=c("native_vs_cultivated")
  } else {
    pair=c("salmiana_vs_tequilana","salmiana_vs_cacti", "cacti_vs_tequilana")
  }
  
  #select elements of the comparison
  sub.ova=data.frame(row.names = sort(unique(path.data$pathway)))
  for (pairs in pair){
    
    ova=read.delim(file = paste(main.dir,comparison,pairs,"OVA.txt",sep = "/"))
    ova=ova[sort(rownames(ova)),]
    colnames(ova)=paste(comparison,pairs,colnames(ova),sep = "-")
    sub.ova=cbind(sub.ova,ova)
    
    
  }
 ova.result= cbind(ova.result,sub.ova)
}

ova.result[ova.result>0.05]=1
ova.result[ova.result<=0.05]=0
ova.result=ova.result[,grep("ALL",colnames(ova.result),value = T, invert = T)]
ova.result=ova.result[,grep("SOI",colnames(ova.result),value = T, invert = T)]
ova.result_2=ova.result[rowSums(ova.result)!=18,]

write.table(ova.result_2,file = "ova_summary.txt", sep = "\t", col.names = T, row.names = T)

## PLOT SUMMARY ##

library(pheatmap)
pal=colorRampPalette(c("black","grey80"))

anot=data.frame(row.names = colnames(data),
                Overrepresented_in=c(rep(c("A.salmiana","A.tequilana","A.salmiana","Cacti","Cacti","A.tequilana"),2),
                                     rep(c("Phyllosphere","Rhizosphere","Phyllosphere","Soil","Rhizosphere","Soil"),2)),
                Comparison_between=c(rep("salmiana_vs_tequilana",2),rep("salmiana_vs_cacti",2),rep("cacti_vs_tequilana",2),
                                     rep("salmiana_vs_tequilana",2),rep("salmiana_vs_cacti",2),rep("cacti_vs_tequilana",2),
                                     rep("phyllosphere_vs_rhizosphere",2),rep("phyllosphere_vs_soil",2),rep("rhizosphere_vs_soil",2),
                                     rep("phyllosphere_vs_rhizosphere",2),rep("phyllosphere_vs_soil",2),rep("rhizosphere_vs_soil",2)),
                Analysis_with=c(rep("Phyllosphere",6),rep("Rhizosphere",6),
                                rep("Native",6),rep("Cultivated",6)))

anot_col=list(Overrepresented_in=c(A.salmiana="darkorange4",
                                   A.tequilana="darkorange3",
                                   Cacti="darkorange1",
                                   Phyllosphere="yellow4",
                                   Rhizosphere="yellow3",
                                   Soil="yellow1"),
              Comparison_between=c(salmiana_vs_tequilana="springgreen4",
                                   salmiana_vs_cacti="springgreen3",
                                   cacti_vs_tequilana="springgreen1",
                                   phyllosphere_vs_rhizosphere="steelblue4",
                                   phyllosphere_vs_soil="steelblue3",
                                   rhizosphere_vs_soil="steelblue1"),
              Analysis_with=c(Phyllosphere="hotpink3",
                              Rhizosphere="hotpink1",
                              Native="plum3",
                              Cultivated="plum1")
              
)


anot=data.frame(row.names = colnames(ova.result_2),
                Overrepresented_in=c(rep(c("A.salmiana","A.tequilana","A.salmiana","Cacti","Cacti","A.tequilana"),1),
                                     rep(c("Phyllosphere","Rhizosphere","Phyllosphere","Soil","Rhizosphere","Soil"),2)),
                Comparison_between=c(rep("salmiana_vs_tequilana",2),rep("salmiana_vs_cacti",2),rep("cacti_vs_tequilana",2),
                                     rep("phyllosphere_vs_rhizosphere",2),rep("phyllosphere_vs_soil",2),rep("rhizosphere_vs_soil",2),
                                     rep("phyllosphere_vs_rhizosphere",2),rep("phyllosphere_vs_soil",2),rep("rhizosphere_vs_soil",2)),
                Analysis_with=c(rep("Phyllosphere",6),
                                rep("Native",6),rep("Cultivated",6)))

anot_col=list(Overrepresented_in=c(A.salmiana="darkorange4",
                                   A.tequilana="darkorange3",
                                   Cacti="darkorange1",
                                   Phyllosphere="yellow4",
                                   Rhizosphere="yellow3",
                                   Soil="yellow1"),
              Comparison_between=c(salmiana_vs_tequilana="springgreen4",
                                   salmiana_vs_cacti="springgreen3",
                                   cacti_vs_tequilana="springgreen1",
                                   phyllosphere_vs_rhizosphere="steelblue4",
                                   phyllosphere_vs_soil="steelblue3",
                                   rhizosphere_vs_soil="steelblue1"),
              Analysis_with=c(Phyllosphere="hotpink3",
                              Native="plum3",
                              Cultivated="plum1")
              
)

# SAVE plot
png(filename = paste(main.dir,"ova.summary_tutorial.png",sep = "/"), width = 320, height = 250, units = "mm", pointsize = 15, bg = "white", res=300, type="cairo")
pheatmap(ova.result_2,
         color = pal(2),
         border_color = "white", 
         cluster_rows = T, cluster_cols = T,
         cellwidth = 8, cellheight =8, fontsize = 9,
         annotation_col = anot,
         annotation_colors = anot_col,
         cutree_cols = 3, cutree_rows = 8,
         labels_col = NA,
         clustering_method = "average")
dev.off()



