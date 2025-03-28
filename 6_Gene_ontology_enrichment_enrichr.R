####EnrichR####
setwd("G:/Headstone/Enh_priming/")
library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr")   
}
library(ggplot2)

#Listing all possible databases
if(websiteLive) dbs = listEnrichrDbs()

#function for adding missing terms when plotting nonprim vs prim groupings
add.terms=function(in.df,terms=tog.terms){
  sub.df=in.df[in.df$Term%in%terms,]
  not.in.terms=terms[!terms%in%sub.df$Term]  
  if(length(not.in.terms)>0){
    not.in.df=as.data.frame(matrix(rep(0),nrow = length(not.in.terms),ncol = ncol(sub.df)))
    colnames(not.in.df)=colnames(sub.df)
    not.in.df$Term=not.in.terms
    not.in.df$Adjusted.P.value=NA
    return(rbind(sub.df,not.in.df))}else{
      return(sub.df)
    }
}


dbs=c("GO_Biological_Process_2023",
  "SynGO_2022",
  "WikiPathway_2021_Human",
  "ENCODE_Histone_Modifications_2015",
  "KEGG_2021_Human")


lookup=read.delim("./Annotations/Human/ENSEMBL_IDs_to_gene_names.txt")
colnames(lookup)=c("Gene.ID","Gene.Name","NCBI_description")


enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(enh.lookup)=c("Enh.ID","Enh.class")
classes=unique(enh.lookup$Enh.class)
prox=read.delim("./Annotations/Human/Prox_genes/Nearest_expressed_gene_20kb.txt",h=F)
colnames(prox)=c("Enh.ID","Gene.ID")
hic=read.delim("./Annotations/Human/HiC_genes/hESC_HiC_genes.txt",h=F)
colnames(hic)=c("Enh.ID","Gene.ID")
comb=rbind(prox,hic)
comb=merge(comb,enh.lookup,by="Enh.ID")
comb=merge(comb,lookup,by="Gene.ID")
head(comb)


coi="NPC_nonprim"
for(n in classes){
coi=n
genes=unique(comb$Gene.Name[comb$Enh.class==coi])
set.size=length(genes)
enriched=enrichr(genes,dbs)
pdf(paste0("./Plots/Gene_ontologies/Human_",coi,"_genes.pdf"))
for(i in 1:length(dbs)){
df=enriched[[i]]
df$num_genes=as.numeric(unlist(lapply(df$Overlap,function(x){unlist(strsplit(x,"/"))[1]})))
df$gene.ratio=df$num_genes/set.size
if(min(df$Adjusted.P.value)>0.05){
df=df[order(df$num_genes,decreasing = T),]
}
df$Term=factor(df$Term,levels=rev(df$Term))


print(ggplot(df[1:5,])+geom_bar(aes(y=Term,weight=gene.ratio,fill=Adjusted.P.value))+
  scale_fill_gradient(low = "red",high="blue",limits=c(0.001,0.05))+xlab("Gene ratio")+ggtitle(paste(coi,"associated genes : ",dbs[i]))+theme_bw())
hea}
dev.off()
}

#ggplot(df,aes(x))


prim.classes=classes[grepl("_prim",classes)]
for(n in prim.classes){
prim.coi=n
non.coi=gsub("_prim","_nonprim",prim.coi)
coi=gsub("_prim","",prim.coi)
prim.genes=unique(comb$Gene.Name[comb$Enh.class==prim.coi])
prim.set.size=length(prim.genes)
prim.enriched=enrichr(prim.genes,dbs)
non.genes=unique(comb$Gene.Name[comb$Enh.class==non.coi])
non.set.size=length(non.genes)
non.enriched=enrichr(non.genes,dbs)

pdf(paste0("./Plots/Gene_ontologies/",coi,"_comparisons.pdf"))
for(i in 1:length(dbs)){
df=prim.enriched[[i]]
df$num_genes=as.numeric(unlist(lapply(df$Overlap,function(x){unlist(strsplit(x,"/"))[1]})))
df$gene.ratio=df$num_genes/prim.set.size
if(min(df$Adjusted.P.value)>0.05){
  df=df[order(df$num_genes,decreasing = T),]
}
df$Term=factor(df$Term,levels=rev(df$Term))
prim.terms=df$Term[1:5]


non.df=non.enriched[[i]]
non.df$num_genes=as.numeric(unlist(lapply(non.df$Overlap,function(x){unlist(strsplit(x,"/"))[1]})))
non.df$gene.ratio=non.df$num_genes/non.set.size
if(min(non.df$Adjusted.P.value)>0.05){
  non.df=non.df[order(non.df$num_genes,decreasing = T),]
}
non.df$Term=factor(non.df$Term,levels=rev(non.df$Term))
non.terms=non.df$Term[1:5]

tog.terms=as.character(c(prim.terms,non.terms))

add.terms=function(in.df,terms=tog.terms){
  sub.df=in.df[in.df$Term%in%terms,]
  not.in.terms=terms[!terms%in%sub.df$Term]  
  if(length(not.in.terms)>0){
  not.in.df=as.data.frame(matrix(rep(0),nrow = length(not.in.terms),ncol = ncol(sub.df)))
  colnames(not.in.df)=colnames(sub.df)
  not.in.df$Term=not.in.terms
  not.in.df$Adjusted.P.value=NA
  return(rbind(sub.df,not.in.df))}else{
    return(sub.df)
  }
  }
prim.sub=add.terms(df)
prim.sub$Term=factor(prim.sub$Term,levels=rev(prim.sub$Term))
non.sub=add.terms(non.df)
non.sub$Term=factor(non.sub$Term,levels=rev(prim.sub$Term))

col.scale=scale_fill_gradient(low = "red",high="blue",limits=c(0,ceiling((max(c(prim.sub$gene.ratio,non.sub$gene.ratio))*10))/10))

print(ggplot(prim.sub)+geom_bar(aes(y=Term,weight=-log10(Adjusted.P.value),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+scale_y_discrete(labels=label_wrap(50))+
  xlab("Adjusted P Value (-log10(p))")+ggtitle(paste(prim.coi,"associated genes : ",dbs[i]))+theme_bw()+col.scale)
print(ggplot(non.sub)+geom_bar(aes(y=Term,weight=-log10(Adjusted.P.value),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+scale_y_discrete(labels=label_wrap(50))+
  xlab("Adjusted P Value (-log10(p))")+ggtitle(paste(non.coi,"associated genes : ",dbs[i]))+theme_bw()+col.scale)

}
dev.off()
}



####Figure plots####



prim.coi="ME_prim"
chosen.terms=c("Adherens junction",
               "(GO:0007411)",
               "(GO:0007399)",
               "(GO:0007179)",
               "(GO:0071902)",
               "Rap1 signaling pathway",
               "(GO:0007409)",
               "Endoderm differentiation WP2853")
get.specific.plots=function(prim.coi,chosen.terms){
non.coi=gsub("_prim","_nonprim",prim.coi)
coi=gsub("_prim","",prim.coi)
pdf(paste0("./Plots/Gene_ontologies/Human/",coi,"_figure_plot.pdf"))

prim.genes=unique(comb$Gene.Name[comb$Enh.class==prim.coi])
prim.set.size=length(prim.genes)
prim.enriched=enrichr(prim.genes,dbs)
non.genes=unique(comb$Gene.Name[comb$Enh.class==non.coi])
non.set.size=length(non.genes)
non.enriched=enrichr(non.genes,dbs)


df.rows=list()
for(i in 1:length(dbs)){
  df=prim.enriched[[i]]
  df$num_genes=as.numeric(unlist(lapply(df$Overlap,function(x){unlist(strsplit(x,"/"))[1]})))
  df$gene.ratio=df$num_genes/prim.set.size
  df.rows[[i]]=df}
tog.df=do.call(rbind,df.rows)
tog.df=tog.df[grepl(pattern = paste(chosen.terms,collapse = "|"),tog.df$Term),]


df.rows=list()
for(i in 1:length(dbs)){
  df=non.enriched[[i]]
  df$num_genes=as.numeric(unlist(lapply(df$Overlap,function(x){unlist(strsplit(x,"/"))[1]})))
  df$gene.ratio=df$num_genes/non.set.size
  df.rows[[i]]=df}
tog.df2=do.call(rbind,df.rows)
tog.df2=tog.df2[grepl(pattern = paste(chosen.terms,collapse = "|"),tog.df2$Term),]

tog.df=tog.df[order(tog.df$Adjusted.P.value,decreasing = T),]
tog.df$Term=factor(tog.df$Term,levels=tog.df$Term)

tog.df2$Term=factor(tog.df2$Term,levels=tog.df$Term)

col.scale=scale_fill_gradient(low = "red",high="blue",limits=c(0,ceiling((max(c(prim.sub$gene.ratio,non.sub$gene.ratio))*10))/10))
print(ggplot(tog.df)+geom_bar(aes(y=Term,weight=-log10(Adjusted.P.value),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+
        xlab("Adjusted P Value (-log10(p))")+ggtitle(paste(prim.coi,"associated genes"))+theme_bw()+col.scale)
print(ggplot(tog.df2)+geom_bar(aes(y=Term,weight=-log10(Adjusted.P.value),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+
        xlab("Adjusted P Value (-log10(p))")+ggtitle(paste(non.coi,"associated genes"))+theme_bw()+col.scale)
dev.off()
}
get.specific.plots("ME_prim",c("Adherens junction",
                               "(GO:0007411)",
                               "(GO:0007399)",
                               "(GO:0007179)",
                               "(GO:0071902)",
                               "Rap1 signaling pathway",
                               "(GO:0007409)",
                               "Endoderm differentiation WP2853"))
get.specific.plots("NPC_prim",c("(GO:0007399)",
                                "(GO:0007409)",
                                "(GO:0099092)",
                                "(GO:0098962)",
                                "WP2858",
                                "WP4249",
                                "WP2853",
                                "Hedgehog signaling pathway"))
get.specific.plots("MSC_prim","04211")
get.specific.plots("7W_brain_prim",c("GO:0007399","GO:0060828","WP2064","Adherens junction","Glycosaminoglycan biosynthesis","Parathyroid hormone synthesis","Rap1 signaling pathway","Oxytocin signaling pathway","MAPK signaling pathway","Neurotrophin signaling pathway"))



####Mouse enrichments####


setwd("G:/Headstone/Enh_priming/")
library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
library(ggplot2)
#install.packages("stringr")
#library(stringr)

#"GO_Biological_Process_2023"
#"SynGO_2022"
#"WikiPathways_2019_Mouse"
#"KEGG_2019_Mouse"

if(websiteLive) dbs = listEnrichrDbs()
dbs$libraryName[order(dbs$libraryName)]

dbs=c(
  "GO_Biological_Process_2023",
  "SynGO_2022",
  "WikiPathways_2019_Mouse",
  "KEGG_2019_Mouse"
)

lookup=read.delim("./Annotations/Mouse/ENSEMBL_IDs_to_gene_names.txt")
colnames(lookup)=c("Gene.ID","Gene.Name","NCBI_ID")

enh.lookup=read.delim("./Annotations/Mouse/Comb_Enh_class_lookup.txt",h=F)
colnames(enh.lookup)=c("Enh.ID","Enh.class")
classes=unique(enh.lookup$Enh.class)

prox=read.delim("./Annotations/Mouse/Prox_genes/Nearest_expressed_gene_20kb.txt",h=F)
colnames(prox)=c("Enh.ID","Gene.ID")
hic=read.delim("./Annotations/Mouse/HiC_genes/EpiLC_HiC_genes.txt",h=F)
colnames(hic)=c("Enh.ID","Gene.ID")
hic2=read.delim("./Annotations/Mouse/HiC_genes/ESC_HiC_genes.txt",h=F)
colnames(hic2)=c("Enh.ID","Gene.ID")
hic=rbind(hic,hic2)

#hic.comb=merge(hic,enh.lookup,by="Enh.ID")
#hic.comb=merge(hic.comb,lookup,by="Gene.ID")
#prox.comb=merge(prox,enh.lookup,by="Enh.ID")
#prox.comb=merge(prox.comb,lookup,by="Gene.ID")

comb=rbind(prox,hic)
comb=merge(comb,enh.lookup,by="Enh.ID")
comb=merge(comb,lookup,by="Gene.ID")

library(scales)
classes

coi="mECT_prim"
classes=classes[grepl("prim",classes)]
for(n in classes){
  coi=n
  genes=unique(comb$Gene.Name[comb$Enh.class==coi])
  set.size=length(genes)
  enriched=enrichr(genes,dbs)
  
  pdf(paste0("./Plots/Gene_ontologies/Mouse/",coi,"_genes.pdf"))
  for(i in 1:length(dbs)){
    df=enriched[[i]]
    df$num_genes=as.numeric(unlist(lapply(df$Overlap,function(x){unlist(strsplit(x,"/"))[1]})))
    df$gene.ratio=df$num_genes/set.size
    
   # if(min(df$Adjusted.P.value)>0.05){
    # df=df[order(df$num_genes,decreasing = T),]
   # }
    df$Term=factor(df$Term,levels=rev(df$Term))
    #df.sub=df[df$Adjusted.P.value<0.2,]
    df.sub=df[1:10,]
    #print(ggplot(df.sub)+geom_bar(aes(y=Term,weight=gene.ratio,fill=Adjusted.P.value))+
     #       scale_fill_gradient(low = "red",high="blue",limits=c(0,0.05))+xlab("Gene ratio")+
    #        scale_y_discrete(labels=label_wrap(50))+
    #       geom_text(aes(y=Term, x=gene.ratio,label=Adjusted.P.value),hjust=0)+
     #       ggtitle(paste(coi,"associated genes : ",dbs[i]))+theme_bw())
    
    col.scale=scale_fill_gradient(low = "red",high="blue",limits=c(0,ceiling((max(c(df.sub$gene.ratio))*10))/10))
    
    print(ggplot(df.sub)+geom_bar(aes(y=Term,weight=-log10(Adjusted.P.value),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+
            scale_y_discrete(labels=label_wrap(50))+xlab("Adjusted P Value (-log10(p))")+
            ggtitle(paste(coi,"associated genes : ",dbs[i]))+theme_bw()+col.scale)
    
  }
  dev.off()
}

coi="heart_prim"
classes=classes[grepl("prim",classes)]
for(n in classes){
  coi=n
  genes=unique(prox.comb$Gene.Name[prox.comb$Enh.class==coi])
  set.size=length(genes)
  enriched=enrichr(genes,dbs)
  
  pdf(paste0("./Plots/Gene_ontologies/Mouse/Prox/",coi,"_genes.pdf"))
  for(i in 1:length(dbs)){
    df=enriched[[i]]
    df$num_genes=as.numeric(unlist(lapply(df$Overlap,function(x){unlist(strsplit(x,"/"))[1]})))
    df$gene.ratio=df$num_genes/set.size
    
    # if(min(df$Adjusted.P.value)>0.05){
    # df=df[order(df$num_genes,decreasing = T),]
    # }
    df$Term=factor(df$Term,levels=rev(df$Term))
    
    
    print(ggplot(df[1:5,])+geom_bar(aes(y=Term,weight=gene.ratio,fill=Adjusted.P.value))+
            scale_fill_gradient(low = "red",high="blue",limits=c(0,0.05))+xlab("Gene ratio")+
            scale_y_discrete(labels=label_wrap(50))+
            geom_text(aes(y=Term, x=gene.ratio,label=Adjusted.P.value),hjust=0)+
            ggtitle(paste(coi,"associated genes : ",dbs[i]))+theme_bw())
  }
  dev.off()
}


library(clusterProfiler)
library(org.Mm.eg.db)

coi="heart_nonprim"
for(n in classes){
  coi=n
  genes=unique(comb$Gene.ID[comb$Enh.class==coi])
  set.size=length(genes)
  GO=enrichGO(genes,OrgDb=org.Mm.eg.db,keyType="ENSEMBL",ont="BP",pAdjustMethod = "BH")
  df=as.data.frame(GO)[,c(1:7,9)]
  if(nrow(df>0)){
  df$gene.ratio=df$Count/set.size
  df$Term=paste0(df$Description," : ",df$ID)
  df$Term=factor(df$Term,levels=rev(df$Term))
  #pdf(paste0("./Plots/Gene_ontologies/Mouse/ClusterProfiler/",coi,"_genes.pdf"))
  print(ggplot(df)+geom_bar(aes(y=Term,weight=gene.ratio,fill=p.adjust))+
          scale_fill_gradient(low = "red",high="blue",limits=c(0,0.05))+xlab("Gene ratio")+ggtitle(paste(coi,"associated genes : GO"))+theme_bw())
  #dev.off()
  }
  }


embryonic.genes=comb$Gene.ID[comb$Enh.class%in%c("mESC_enh","mEpiLC_specific")]
embryonic.genes.names=comb$Gene.Name[comb$Enh.class%in%c("mESC_enh","mEpiLC_specific")]

coi="brain_prim"
genes=unique(comb$Gene.ID[comb$Enh.class%in%c("forebrain_prim","midbrain_prim","hindbrain_prim")])
genes.sub=genes[!genes%in%embryonic.genes]
set.size=length(genes.sub)
GO=enrichGO(genes.sub,OrgDb=org.Mm.eg.db,keyType="ENSEMBL",ont="BP",pAdjustMethod = "BH")
df=as.data.frame(GO)[,c(1:7,9)]
df$gene.ratio=df$Count/set.size
df$Term=paste0(df$Description," : ",df$ID)
df$Term=factor(df$Term,levels=rev(df$Term))
print(ggplot(df)+geom_bar(aes(y=Term,weight=gene.ratio,fill=p.adjust))+
        scale_fill_gradient(low = "red",high="blue",limits=c(0,0.05))+xlab("Gene ratio")+ggtitle(paste(coi,"associated genes : GO"))+theme_bw())

genes=unique(comb$Gene.Name[comb$Enh.class%in%c("forebrain_nonprim","midbrain_nonprim","hindbrain_nonprim")])
genes.sub=genes
genes.sub=genes[!genes%in%embryonic.genes.names]
set.size=length(genes.sub)

enriched=enrichr(genes.sub,dbs)
for(i in 1:length(dbs)){
  df=enriched[[i]]
  df$num_genes=as.numeric(unlist(lapply(df$Overlap,function(x){unlist(strsplit(x,"/"))[1]})))
  df$gene.ratio=df$num_genes/set.size
  
  # if(min(df$Adjusted.P.value)>0.05){
  # df=df[order(df$num_genes,decreasing = T),]
  # }
  df$Term=factor(df$Term,levels=rev(df$Term))
  #df.sub=df[df$Adjusted.P.value<0.2,]
  df.sub=df[1:10,]
  #print(ggplot(df.sub)+geom_bar(aes(y=Term,weight=gene.ratio,fill=Adjusted.P.value))+
  #       scale_fill_gradient(low = "red",high="blue",limits=c(0,0.05))+xlab("Gene ratio")+
  #        scale_y_discrete(labels=label_wrap(50))+
  #       geom_text(aes(y=Term, x=gene.ratio,label=Adjusted.P.value),hjust=0)+
  #       ggtitle(paste(coi,"associated genes : ",dbs[i]))+theme_bw())
  
  col.scale=scale_fill_gradient(low = "red",high="blue",limits=c(0,ceiling((max(c(df.sub$gene.ratio))*10))/10))
  
  print(ggplot(df.sub)+geom_bar(aes(y=Term,weight=-log10(Adjusted.P.value),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+
          scale_y_discrete(labels=label_wrap(50))+xlab("Adjusted P Value (-log10(p))")+
          ggtitle(paste(coi,"associated genes : ",dbs[i]))+theme_bw()+col.scale)
  
}



####Analysis of Mouse gene sets using http://geneontology.org/ ####

classes=classes[grepl("prim",classes)]
for(i in classes){
genes=unique(comb$Gene.Name[comb$Enh.class%in%c(i)])
write.table(genes,paste0("./Data/mouse/GO_enrichment/Gene_names/",i,"_gene_names.txt"),sep="\t",col.names = F,row.names = F,quote = F)
}

#forebrain GO plots
coi="forebrain"
GO.prim=read.delim("./Data/mouse/GO_enrichment/forebrain_prim.txt",skip =11)
num.prim=as.numeric(gsub("\\.","",gsub("upload_1..","",colnames(GO.prim)[3])))
GO.nonprim=read.delim("./Data/mouse/GO_enrichment/forebrain_nonprim.txt",skip=11)
num.nonprim=as.numeric(gsub("\\.","",gsub("upload_1..","",colnames(GO.nonprim)[3])))

processes=unique(GO.prim$GO.biological.process.complete)
GO.nonprim=GO.nonprim[GO.nonprim$GO.biological.process.complete%in%processes,]
head(GO.nonprim)

GO.prim$gene.ratio=GO.prim[,3]/num.prim
GO.nonprim$gene.ratio=GO.nonprim[,3]/num.nonprim

pdf("./Plots/Gene_ontologies/Mouse/Panther/forebrain.pdf")
col.scale=scale_fill_gradient(low = "red",high="blue",limits=c(0,ceiling((max(c(GO.prim$gene.ratio,GO.nonprim$gene.ratio))*10))/10))
GO.prim=GO.prim[order(GO.prim$upload_1..FDR.,decreasing = F),]
GO.prim$GO.biological.process.complete=factor(GO.prim$GO.biological.process.complete,levels=rev(GO.prim$GO.biological.process.complete))
print(ggplot(GO.prim)+geom_bar(aes(y=GO.biological.process.complete,weight=-log10(upload_1..FDR.),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+
        scale_y_discrete(labels=label_wrap(50))+xlab("Adjusted P Value (-log10(p))")+
        ggtitle(paste0(coi,"_prim associated genes"))+theme_bw()+col.scale)
GO.nonprim$GO.biological.process.complete=factor(GO.nonprim$GO.biological.process.complete,levels=rev(GO.prim$GO.biological.process.complete))
print(ggplot(GO.nonprim)+geom_bar(aes(y=GO.biological.process.complete,weight=-log10(upload_1..FDR.),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+
        scale_y_discrete(labels=label_wrap(50))+xlab("Adjusted P Value (-log10(p))")+
        ggtitle(paste0(coi,"_nonprim associated genes"))+theme_bw()+col.scale)
dev.off()

#mECT plots
coi="mECT"
GO.prim=read.delim("./Data/mouse/GO_enrichment/mECT_prim.txt",skip =11)
num.prim=as.numeric(gsub("\\.","",gsub("upload_1..","",colnames(GO.prim)[3])))
GO.nonprim=read.delim("./Data/mouse/GO_enrichment/mECT_nonprim.txt",skip=11)
num.nonprim=as.numeric(gsub("\\.","",gsub("upload_1..","",colnames(GO.nonprim)[3])))

processes=unique(GO.prim$GO.biological.process.complete)
GO.nonprim=GO.nonprim[GO.nonprim$GO.biological.process.complete%in%processes,]
GO.prim$gene.ratio=GO.prim[,3]/num.prim
GO.nonprim$gene.ratio=GO.nonprim[,3]/num.nonprim

pdf("./Plots/Gene_ontologies/Mouse/Panther/mECT.pdf")
col.scale=scale_fill_gradient(low = "red",high="blue",limits=c(0,ceiling((max(c(GO.prim$gene.ratio,GO.nonprim$gene.ratio))*10))/10))
GO.prim=GO.prim[order(GO.prim$upload_1..FDR.,decreasing = F),]
GO.prim$GO.biological.process.complete=factor(GO.prim$GO.biological.process.complete,levels=rev(GO.prim$GO.biological.process.complete))
print(ggplot(GO.prim)+geom_bar(aes(y=GO.biological.process.complete,weight=-log10(upload_1..FDR.),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+
        scale_y_discrete(labels=label_wrap(50))+xlab("Adjusted P Value (-log10(p))")+
        ggtitle(paste0(coi,"_prim associated genes"))+theme_bw()+col.scale)
GO.nonprim$GO.biological.process.complete=factor(GO.nonprim$GO.biological.process.complete,levels=rev(GO.prim$GO.biological.process.complete))
print(ggplot(GO.nonprim)+geom_bar(aes(y=GO.biological.process.complete,weight=-log10(upload_1..FDR.),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+
        scale_y_discrete(labels=label_wrap(50))+xlab("Adjusted P Value (-log10(p))")+
        ggtitle(paste0(coi,"_nonprim associated genes"))+theme_bw()+col.scale)
dev.off()

#mEND plots
coi="mEND"
GO.prim=read.delim("./Data/mouse/GO_enrichment/mEND_prim.txt",skip =11)
num.prim=as.numeric(gsub("\\.","",gsub("upload_1..","",colnames(GO.prim)[3])))
GO.nonprim=read.delim("./Data/mouse/GO_enrichment/mEND_nonprim.txt",skip=11)
num.nonprim=as.numeric(gsub("\\.","",gsub("upload_1..","",colnames(GO.nonprim)[3])))

processes=unique(GO.prim$GO.biological.process.complete)
GO.nonprim=GO.nonprim[GO.nonprim$GO.biological.process.complete%in%processes,]
GO.prim$gene.ratio=GO.prim[,3]/num.prim
GO.nonprim$gene.ratio=GO.nonprim[,3]/num.nonprim

pdf("./Plots/Gene_ontologies/Mouse/Panther/mEND.pdf")
col.scale=scale_fill_gradient(low = "red",high="blue",limits=c(0,ceiling((max(c(GO.prim$gene.ratio,GO.nonprim$gene.ratio))*10))/10))
GO.prim=GO.prim[order(GO.prim$upload_1..FDR.,decreasing = F),]
GO.prim$GO.biological.process.complete=factor(GO.prim$GO.biological.process.complete,levels=rev(GO.prim$GO.biological.process.complete))
print(ggplot(GO.prim)+geom_bar(aes(y=GO.biological.process.complete,weight=-log10(upload_1..FDR.),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+
        scale_y_discrete(labels=label_wrap(50))+xlab("Adjusted P Value (-log10(p))")+
        ggtitle(paste0(coi,"_prim associated genes"))+theme_bw()+col.scale)
GO.nonprim$GO.biological.process.complete=factor(GO.nonprim$GO.biological.process.complete,levels=rev(GO.prim$GO.biological.process.complete))
print(ggplot(GO.nonprim)+geom_bar(aes(y=GO.biological.process.complete,weight=-log10(upload_1..FDR.),fill=gene.ratio))+xlim(c(0,5))+geom_vline(xintercept = -log10(0.05))+
        scale_y_discrete(labels=label_wrap(50))+xlab("Adjusted P Value (-log10(p))")+
        ggtitle(paste0(coi,"_nonprim associated genes"))+theme_bw()+col.scale)
dev.off()




####Getting gene ID files for output####
setwd("G:/Headstone/Enh_priming/")

lookup=read.delim("./Annotations/Mouse/ENSEMBL_IDs_to_gene_names.txt")
colnames(lookup)=c("Gene.ID","Gene.Name","NCBI_ID")
enh.lookup=read.delim("./Annotations/Mouse/Comb_Enh_class_lookup.txt",h=F)
colnames(enh.lookup)=c("Enh.ID","Enh.class")
classes=unique(enh.lookup$Enh.class)

prox=read.delim("./Annotations/Mouse/Prox_genes/Nearest_expressed_gene_20kb.txt",h=F)
colnames(prox)=c("Enh.ID","Gene.ID")
hic=read.delim("./Annotations/Mouse/HiC_genes/EpiLC_HiC_genes.txt",h=F)
colnames(hic)=c("Enh.ID","Gene.ID")
hic2=read.delim("./Annotations/Mouse/HiC_genes/ESC_HiC_genes.txt",h=F)
colnames(hic2)=c("Enh.ID","Gene.ID")
hic=rbind(hic,hic2)


comb=rbind(prox,hic)
comb=merge(comb,enh.lookup,by="Enh.ID")
comb=merge(comb,lookup,by="Gene.ID")


classes=classes[grepl("prim",classes)]
for(n in classes){
  coi=n
  genes=unique(comb$Gene.Name[comb$Enh.class==coi])
  write.table(genes,paste0("./Annotations/Mouse/Associated_genes/",n,"_genes.txt"),sep="\t",col.names = F,row.names = F,quote = F)
}




