####Human####

setwd("F:/Enh_priming/")

library(ggplot2)
library(tidyverse)
library(qlcMatrix)
library(ggpubr)


{
##Loading protein coding genes

tss=read.delim("./Annotations/Human/ENSEMBL_coding_gene_TSS.txt",h=F)
gene.ids=unique(tss[,4])

rel.classes=c("NPC_prim"  , "NPC_nonprim","ME_prim" ,"ME_nonprim","MSC_prim", "MSC_nonprim")

##Function to add coding genes which were not covered by a dataset and sets them to the minimal observed value
add.not.in=function(df){
  not.in.ids=gene.ids[!gene.ids%in%df[,1]]
  min.val=min(df[,2:ncol(df)]) 
  not.in.mat=as.data.frame(matrix(data=rep(min.val),ncol = ncol(df),nrow=length(not.in.ids)))
  not.in.mat[,1]=not.in.ids
  colnames(not.in.mat)=colnames(df)
  return(rbind(df,not.in.mat))
}


##Loading and formating all relevant expression datasets
xie=read.delim("./Data/human/RNA/Xie_et_al_RNA_logRPKM_all_types.txt")[,c(7,13:16)]
#xie=add.not.in(xie)
xie=pivot_longer(xie,cols = -1, names_to = "Tissue",values_to = "Expr")
colnames(xie)[1]="Gene.ID"
xie.sub=xie[xie$Gene.ID%in%gene.ids,]


####Getting the human enhancer groups and gene associations####
enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(enh.lookup)=c("Enh.ID","Enh.class")
classes=unique(enh.lookup$Enh.class)

prox=read.delim("./Annotations/Human/Prox_genes/Nearest_expressed_gene_20kb.txt",h=F)
colnames(prox)=c("Enh.ID","Gene.ID")
hic=read.delim("./Annotations/Human/HiC_genes/hESC_HiC_genes.txt",h=F)
colnames(hic)=c("Enh.ID","Gene.ID")
comb=rbind(prox,hic)

comb=merge(comb,enh.lookup,by="Enh.ID")


embryonic.genes=unique(comb$Gene.ID[grepl("hESC_specific",comb$Enh.class)])

prescreen=comb[!comb$Gene.ID%in%embryonic.genes,]
prescreen=prescreen[prescreen$Enh.class%in%rel.classes,]


prim.classes=rel.classes[grepl("_prim",rel.classes)]
df.rows=list()
for(i in rel.classes){
  if(grepl("_prim",i)){
    gene.list=unique(comb$Gene.ID[comb$Enh.class==i])
    other.prim.classes=prim.classes[prim.classes!=i]
    other.prim=comb$Gene.ID[comb$Enh.class%in%other.prim.classes]
    df.rows[[paste0(i,".embryonic")]]=data.frame("Enh.class"=rep(paste0(i,".embryonic")),"Gene.ID"=gene.list[gene.list%in%embryonic.genes&!gene.list%in%other.prim])
    gene.list=gene.list[!gene.list%in%embryonic.genes]
    other.prim.classes=prim.classes[prim.classes!=i]
    other.prim=comb$Gene.ID[comb$Enh.class%in%other.prim.classes]
    df.rows[[paste0(i,".specific")]]=data.frame("Enh.class"=rep(i),"Gene.ID"=gene.list[!gene.list%in%other.prim])
  }else   if(grepl("_nonprim",i)){
    gene.list=unique(comb$Gene.ID[comb$Enh.class==i])
    #gene.list=gene.list[!gene.list%in%embryonic.genes]
    prim.list=comb$Gene.ID[grepl("_prim",comb$Enh.class)]
    gene.list=gene.list[!gene.list%in%prim.list]
    df.rows[[i]]=data.frame("Enh.class"=rep(i),"Gene.ID"=gene.list)
  }
}

df.tog=do.call(rbind,df.rows)
}

head(enh.comb)

enh.comb=merge(xie,df.tog,by="Gene.ID")
pdf("./Expression_revisit/all_genes_screened_associations_logRPM.pdf")
ggplot(enh.comb[enh.comb$Enh.class=="NPC_prim",],aes(x=Tissue,y=Expr))+geom_violin()+geom_boxplot(width=0.1)+theme_bw()+ggtitle("NPC_prim")+ylim(c(-13,10))
ggplot(enh.comb[enh.comb$Enh.class=="NPC_nonprim",],aes(x=Tissue,y=Expr))+geom_violin()+geom_boxplot(width=0.1)+theme_bw()+ggtitle("NPC_nonprim")+ylim(c(-13,10))
ggplot(enh.comb[enh.comb$Enh.class=="ME_prim",],aes(x=Tissue,y=Expr))+geom_violin()+geom_boxplot(width=0.1)+theme_bw()+ggtitle("ME_prim")+ylim(c(-13,10))
ggplot(enh.comb[enh.comb$Enh.class=="ME_nonprim",],aes(x=Tissue,y=Expr))+geom_violin()+geom_boxplot(width=0.1)+theme_bw()+ggtitle("ME_nonprim")+ylim(c(-13,10))
ggplot(enh.comb[enh.comb$Enh.class=="MSC_prim",],aes(x=Tissue,y=Expr))+geom_violin()+geom_boxplot(width=0.1)+theme_bw()+ggtitle("MSC_prim")+ylim(c(-13,10))
ggplot(enh.comb[enh.comb$Enh.class=="MSC_nonprim",],aes(x=Tissue,y=Expr))+geom_violin()+geom_boxplot(width=0.1)+theme_bw()+ggtitle("MSC_nonprim")+ylim(c(-13,10))
dev.off()

t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="NPC_prim"&enh.comb$Tissue%in%c("H1","NPC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="NPC_prim"&enh.comb$Tissue%in%c("H1","MSC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="NPC_prim"&enh.comb$Tissue%in%c("H1","ME"),])

t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="NPC_nonprim"&enh.comb$Tissue%in%c("H1","NPC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="NPC_nonprim"&enh.comb$Tissue%in%c("H1","MSC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="NPC_nonprim"&enh.comb$Tissue%in%c("H1","ME"),])

t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="ME_prim"&enh.comb$Tissue%in%c("H1","NPC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="ME_prim"&enh.comb$Tissue%in%c("H1","MSC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="ME_prim"&enh.comb$Tissue%in%c("H1","ME"),])

t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="ME_nonprim"&enh.comb$Tissue%in%c("H1","NPC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="ME_nonprim"&enh.comb$Tissue%in%c("H1","MSC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="ME_nonprim"&enh.comb$Tissue%in%c("H1","ME"),])


t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="MSC_prim"&enh.comb$Tissue%in%c("H1","NPC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="MSC_prim"&enh.comb$Tissue%in%c("H1","MSC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="MSC_prim"&enh.comb$Tissue%in%c("H1","ME"),])

t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="MSC_nonprim"&enh.comb$Tissue%in%c("H1","NPC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="MSC_nonprim"&enh.comb$Tissue%in%c("H1","MSC"),])
t.test(Expr ~ Tissue, enh.comb[enh.comb$Enh.class=="MSC_nonprim"&enh.comb$Tissue%in%c("H1","ME"),])


setwd("G:/Headstone/Enh_priming/")

library(ggplot2)
library(tidyverse)
library(qlcMatrix)
library(ggpubr)

tss=read.delim("./Annotations/Mouse/mm10_ENSEMBL_coding_gene_TSS.bed",h=F)

{
  ##Load all expression datasets
  esc=read.delim("./Data/mouse/RNA/GSE223575_expression_gene_log2RPKM.txt")[,c(6,8:9,46:47,12:13,52:53)]
  esc$ESC=rowMeans(esc[,2:5])
  esc$EpiLC=rowMeans(esc[,6:9])
  esc=esc[,c(1,10:11)]
  colnames(esc)[1]="Gene.ID"
  
  #Some genes are measured more than once, will take just the highest expression version for each tissue
  esc.esc=esc[,1:2]
  esc.epi=esc[,c(1,3)]
  esc.esc=esc.esc%>%group_by(Gene.ID)%>%slice(which.max(ESC))%>%ungroup()
  esc.epi=esc.epi%>%group_by(Gene.ID)%>%slice(which.max(EpiLC))%>%ungroup()
  esc=merge(esc.esc,esc.epi,by="Gene.ID")
  rm(esc.esc,esc.epi)
  
  #scnmt seq data is "log normalised counts"
  scnmt=read.delim("./Data/mouse/RNA/scNMT_pseudobulk_logcounts.txt")
  scnmt$Gene.ID=rownames(scnmt)
  scnmt=scnmt[,c(8,1:7)]
  
  fetal=read.delim("./Data/mouse/RNA/GSE100685_TPM_matrix.tsv")
  fetal$GeneID=unlist(lapply(fetal$GeneID,function(x){unlist(strsplit(x,"\\."))[1]}))
  fetal$E11_Brain.fore=log(rowMeans(fetal[,c(2,10)]))
  fetal$E11_Brain.mid=log(rowMeans(fetal[,c(48,56)]))
  fetal$E11_Brain.hind=log(rowMeans(fetal[,c(18,26)]))
  fetal$E11_Liver=log(rowMeans(fetal[,c(33,40)]))
  fetal=fetal[,c(1,63:66)]
  colnames(fetal)[1]="Gene.ID"
  
  fetal.heart=read.delim("./Data/mouse/RNA/GSE205165_Processed_data_for_GEO_ADRC_vs._E11.5_CM.txt")
  fetal.heart$E11_Heart=log(rowMeans(fetal.heart[,4:5]))
  fetal.heart$E11_Adipose=log(rowMeans(fetal.heart[,2:3]))
  fetal.heart=fetal.heart[,c(1,6:7)]
  colnames(fetal.heart)[1]="Gene.ID"
  
  gene.ids=unique(c(tss[,4],esc[,1],scnmt[,1],fetal[,1],fetal.heart[,1]))
  
  add.not.measured.genes=function(df){
    not.in.genes=gene.ids[!gene.ids%in%df[,1]]
    not.mat=as.data.frame(matrix(data=rep(-Inf),ncol=ncol(df),nrow=length(not.in.genes)))
    not.mat[,1]=not.in.genes
    colnames(not.mat)=colnames(df)
    new.df=rbind(df,not.mat)
    return(new.df)
  }
  
  esc=add.not.measured.genes(esc)
  scnmt=add.not.measured.genes(scnmt)
  fetal=add.not.measured.genes(fetal)
  fetal.heart=add.not.measured.genes(fetal.heart)
  
  
  somatic=read.delim("./Data/mouse/RNA/ENCODE_somatic_logRPKM.txt")[,c(7,13:17)]
  somatic=add.not.measured.genes(somatic)
  colnames(somatic)[1]="Gene.ID"
  
  df.rows=list()
  for(i in gene.ids){
    df.sub=somatic[somatic$Gene.ID==i,]
    df.rows[[i]]=c(i,apply(as.matrix(df.sub[,-1]),2,max))
  }
  tog=as.data.frame(do.call(rbind,df.rows))
  colnames(tog)=colnames(somatic)
  somatic=tog
  rm(df.rows,tog)
  somatic=somatic[gene.ids,]
  colnames(somatic)[2:6]=paste0("somatic_",colnames(somatic)[2:6])
  for(i in 2:6){
    somatic[,i]=as.numeric(somatic[,i])
  }
  
  rownames(esc)=esc$Gene.ID
  rownames(scnmt)=scnmt$Gene.ID
  rownames(fetal)=fetal$Gene.ID
  rownames(fetal.heart)=fetal.heart$Gene.ID
  esc=esc[gene.ids,]
  scnmt=scnmt[gene.ids,]
  fetal=fetal[gene.ids,]
  fetal.heart=fetal.heart[gene.ids,]
  
  cutoff.esc=1
  cutoff.scnmt=10
  cutoff.fetal=1
  cutoff.fetal.heart=1
  
  expr.esc=as.numeric((rowMax(as.matrix(esc[,2:3]))))>cutoff.esc
  expr.scnmt=as.numeric((rowMax(as.matrix(scnmt[,2:3]))))>cutoff.scnmt
  expr.fetal=as.numeric((rowMax(as.matrix(fetal[,2:3]))))>cutoff.fetal
  expr.fetal.heart=as.numeric((rowMax(as.matrix(fetal.heart[,2:3]))))>cutoff.fetal.heart
  
  expr.genes=gene.ids[expr.esc|expr.scnmt|expr.fetal|expr.fetal.heart]
  rm(expr.esc,expr.scnmt,expr.fetal,expr.fetal.heart,cutoff.esc,cutoff.fetal,cutoff.fetal.heart,cutoff.scnmt)
  
  data=do.call(cbind,list(esc,scnmt[,2:8],fetal[,2:5],fetal.heart[,2:3],somatic[,2:6]))
  data.coding=data[data$Gene.ID%in%tss[,4],]
  data.expr=data[data$Gene.ID%in%expr.genes,]
  
  
  
  data=pivot_longer(data,cols=-1,names_to = "Tissue",values_to="Expression")
  
  data.coding=data[data$Gene.ID%in%tss[,4],]%>%group_by(Tissue)%>%mutate(percentile=percent_rank(Expression)*100)
  data.expr=data[data$Gene.ID%in%expr.genes,]%>%group_by(Tissue)%>%mutate(percentile=percent_rank(Expression)*100)
  data.all=data%>%group_by(Tissue)%>%mutate(percentile=percent_rank(Expression)*100)
  
  ctrl.coding=data.coding
  ctrl.expr=data.expr
  ctrl.all=data.all
  ctrl.coding$Enh.class=ctrl.expr$Enh.class=ctrl.all$Enh.class=rep("Global")
}
{
  ###Loading Enhancer-gene associations
  
  enh.lookup=read.delim("./Annotations/Mouse/Comb_Enh_class_lookup.txt",h=F)
  colnames(enh.lookup)=c("ID","Enh.class")
  classes=unique(enh.lookup$Enh.class)
  
  epi.hic=read.delim("./Annotations/Mouse/HiC_genes/EpiLC_HiC_genes.txt",h=F)
  prox=read.delim("./Annotations/Mouse/Prox_genes/Nearest_expressed_gene_20kb.txt",h=F)
  colnames(epi.hic)=colnames(prox)=c("ID","Gene.ID")
  
  comb=rbind(epi.hic,prox)
  comb=merge(comb,enh.lookup,by="ID")
  
  ###Treating the enhancer-gene associations
  ##Removing Embryonic enhancer associated genes, and genes regulated by more than one lineage of enhancers
  embryonic.genes=unique(comb$Gene.ID[grepl("ESC|EpiLC",comb$Enh.class)])
  
  rel.classes=c("mECT_prim","mECT_nonprim","mEND_prim","mEND_nonprim","mMES_prim","mMES_nonprim","liver_prim","liver_nonprim",    
                "heart_prim","heart_nonprim","forebrain_prim","forebrain_nonprim","midbrain_prim","midbrain_nonprim","hindbrain_prim","hindbrain_nonprim")
  prim.classes=rel.classes[grepl("_prim",rel.classes)]
  df.rows=list()
  for(i in rel.classes){
    if(grepl("_prim",i)){
      gene.list=unique(comb$Gene.ID[comb$Enh.class==i])
      other.prim.classes=prim.classes[prim.classes!=i]
      other.prim=comb$Gene.ID[comb$Enh.class%in%other.prim.classes]
      df.rows[[paste0(i,".embryonic")]]=data.frame("Enh.class"=rep(paste0(i,".embryonic")),"Gene.ID"=gene.list[gene.list%in%embryonic.genes&!gene.list%in%other.prim])
      gene.list=gene.list[!gene.list%in%embryonic.genes]
      other.prim.classes=prim.classes[prim.classes!=i]
      other.prim=comb$Gene.ID[comb$Enh.class%in%other.prim.classes]
      df.rows[[paste0(i,".specific")]]=data.frame("Enh.class"=rep(i),"Gene.ID"=gene.list[!gene.list%in%other.prim])
    }else   if(grepl("_nonprim",i)){
      gene.list=unique(comb$Gene.ID[comb$Enh.class==i])
      # gene.list=gene.list[!gene.list%in%embryonic.genes]
      prim.list=comb$Gene.ID[grepl("_prim",comb$Enh.class)]
      gene.list=gene.list[!gene.list%in%prim.list]
      df.rows[[i]]=data.frame("Enh.class"=rep(i),"Gene.ID"=gene.list)
    }
  }
  df.tog=do.call(rbind,df.rows)
  
  
  non.neural.prim=c("mECT_prim","mEND_prim","mMES_prim","liver_prim","heart_prim")
  brain.prim.genes=unique(comb$Gene.ID[comb$Enh.class=="brain_prim"])
  other.prim.genes=comb$Gene.ID[comb$Enh.class%in%non.neural.prim]
  brain.prim.genes=brain.prim.genes[!brain.prim.genes%in%c(embryonic.genes,other.prim.genes)]
  brain.prim.df=data.frame("Enh.class"=rep("brain_prim"),"Gene.ID"=brain.prim.genes)
  brain.non.genes=unique(comb$Gene.ID[comb$Enh.class=="brain_nonprim"])
  other.prim.genes=comb$Gene.ID[comb$Enh.class%in%prim.classes]
  brain.non.genes=brain.non.genes[!brain.non.genes%in%c(embryonic.genes,other.prim.genes)]
  brain.nonprim.df=data.frame("Enh.class"=rep("brain_nonprim"),"Gene.ID"=brain.non.genes)
  df.tog=do.call(rbind,list(df.tog,brain.prim.df,brain.nonprim.df))
  rm(brain.prim.genes,other.prim.genes,brain.non.genes,brain.nonprim.df,non.neural.prim)
  
  df.emb=data.frame("Enh.class"=rep("Embryonic"),"Gene.ID"=embryonic.genes)
  df.tog=rbind(df.tog,df.emb)
  plot.order=c("ESC","EpiLC","E4.5_Epiblast","E5.5_Epiblast","E6.5_Epiblast","E6.5_Primitive_Streak",
               "E7.5_Ectoderm","E7.5_Endoderm","E7.5_Mesoderm",
               "E11_Brain.fore","E11_Brain.mid","E11_Brain.hind","E11_Liver","E11_Heart","E11_Adipose",
               "somatic_adrenal_gland","somatic_gastrocnemius","somatic_heart","somatic_left_cerebral_cortex","somatic_hippocampus")
  head(enh.coding)
  unique(enh.coding$Tissue)
  ctrl.subset=sample(gene.ids,size=1000,replace=FALSE)
  ctrl.expr=data.expr[data.expr$Gene.ID%in%ctrl.subset,]
  ctrl.expr$Enh.class=rep("Rand.Ctrl")
  
  enh.all=merge(data.all,df.tog,by="Gene.ID")
  enh.expr=merge(data.expr,df.tog,by="Gene.ID")
  enh.coding=merge(data.coding,df.tog,by="Gene.ID")
  
  enh.expr=rbind(enh.expr,ctrl.expr)
  
  enh.all$Tissue=factor(enh.all$Tissue,levels = plot.order)
  enh.expr$Tissue=factor(enh.expr$Tissue,levels = plot.order)
  enh.coding$Tissue=factor(enh.coding$Tissue,levels = plot.order)
  
  
  
}

pdf("./Removing_Embryonic_enhancer_Associated_genes.pdf")

venn=list("ESC/EpiLC Enh"=unique(comb$Gene.ID[grepl("ESC|EpiLC",comb$Enh.class)]),
          "mECT_prim Enh"=unique(comb$Gene.ID[grepl("mECT_prim",comb$Enh.class)]),
          "mEND_prim Enh"=unique(comb$Gene.ID[grepl("mEND_prim",comb$Enh.class)]),
          "mMES_prim Enh"=unique(comb$Gene.ID[grepl("mMES_prim",comb$Enh.class)]))
print(ggvenn(venn))

venn=list("ESC/EpiLC Enh"=unique(comb$Gene.ID[grepl("ESC|EpiLC",comb$Enh.class)]),
          "forebrain_prim Enh"=unique(comb$Gene.ID[grepl("forebrain_prim",comb$Enh.class)]),
          "heart_prim Enh"=unique(comb$Gene.ID[grepl("heart_prim",comb$Enh.class)]),
          "Liver_prim Enh"=unique(comb$Gene.ID[grepl("liver_prim",comb$Enh.class)]))
print(ggvenn(venn))

tissue.plot=c("E4.5_Epiblast","E5.5_Epiblast","E7.5_Ectoderm","E7.5_Endoderm","E7.5_Mesoderm")
to.plot=c("mECT_prim","mECT_prim.embryonic","mEND_prim","mEND_prim.embryonic","mMES_prim","mMES_prim.embryonic")
#ggplot(enh.coding[enh.coding$Enh.class%in%to.plot&enh.coding$Tissue%in%tissue.plot,],aes(x=Enh.class,y=Expression,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw()
print(ggplot(enh.coding[enh.coding$Enh.class%in%to.plot&enh.coding$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())


tissue.plot=c("E11_Brain.fore","E11_Liver","E11_Heart","E11_Adipose")

to.plot=c("forebrain_prim","forebrain_prim.embryonic","heart_prim","heart_prim.embryonic","liver_prim","liver_prim.embryonic")
#ggplot(enh.coding[enh.coding$Enh.class%in%to.plot&enh.coding$Tissue%in%tissue.plot,],aes(x=Enh.class,y=Expression,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw()
print(ggplot(enh.coding[enh.coding$Enh.class%in%to.plot&enh.coding$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())

dev.off()



pdf("./Justification_of_percentile_approach.pdf")
to.plot=c("mMES_prim","mMES_nonprim")
tissue.plot=c("ESC","EpiLC","E4.5_Epiblast","E5.5_Epiblast","E6.5_Epiblast","E6.5_Primitive_Streak",
              "E7.5_Ectoderm","E7.5_Endoderm","E7.5_Mesoderm",
              "E11_Brain.fore","E11_Brain.mid","E11_Brain.hind","E11_Liver","E11_Heart","E11_Adipose")
print(ggplot(enh.coding[enh.coding$Enh.class%in%to.plot&enh.coding$Tissue%in%tissue.plot,],aes(x=Enh.class,y=Expression,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(enh.coding[enh.coding$Enh.class%in%to.plot&enh.coding$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
tissue.plot=c(              "E7.5_Ectoderm","E7.5_Endoderm","E7.5_Mesoderm"              )
print(ggplot(enh.coding[enh.coding$Enh.class%in%to.plot&enh.coding$Tissue%in%tissue.plot,],aes(x=Enh.class,y=Expression,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(enh.coding[enh.coding$Enh.class%in%to.plot&enh.coding$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
dev.off()



##Ploting final comparisons
pdf("./Plots/Gene_expression/Mouse_Prim_nonprim_comparisons_across_dev.pdf")
enh.plot=c("mECT_prim","mEND_prim","mMES_prim","mECT_nonprim","mEND_nonprim","mMES_nonprim")
sub.df=enh.coding[enh.coding$Enh.class%in%enh.plot,]
#sub.df$Enh.class=factor(sub.df$Enh.class,levels=enh.plot)

tissue.plot=c("ESC","EpiLC","E4.5_Epiblast","E5.5_Epiblast","E6.5_Epiblast","E6.5_Primitive_Streak",
              "E7.5_Ectoderm","E7.5_Endoderm","E7.5_Mesoderm",
              "E11_Brain.fore","E11_Brain.mid","E11_Brain.hind","E11_Liver","E11_Heart","E11_Adipose")

print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())

tissue.plot=c("somatic_adrenal_gland","somatic_gastrocnemius","somatic_heart","somatic_left_cerebral_cortex","somatic_hippocampus")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())

enh.plot=c("forebrain_prim","forebrain_nonprim","midbrain_prim","midbrain_nonprim","hindbrain_prim","hindbrain_nonprim")
sub.df=enh.coding[enh.coding$Enh.class%in%enh.plot,]
sub.df$Enh.class=factor(sub.df$Enh.class,levels=enh.plot)

tissue.plot=c("ESC","EpiLC","E4.5_Epiblast","E5.5_Epiblast","E6.5_Epiblast","E6.5_Primitive_Streak",
              "E7.5_Ectoderm","E7.5_Endoderm","E7.5_Mesoderm",
              "E11_Brain.fore","E11_Brain.mid","E11_Brain.hind","E11_Liver","E11_Heart","E11_Adipose")

print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())

tissue.plot=c("somatic_adrenal_gland","somatic_gastrocnemius","somatic_heart","somatic_left_cerebral_cortex","somatic_hippocampus")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())


enh.plot=c("liver_prim","liver_nonprim","heart_prim","heart_nonprim")
sub.df=enh.coding[enh.coding$Enh.class%in%enh.plot,]
sub.df$Enh.class=factor(sub.df$Enh.class,levels=enh.plot)

tissue.plot=c("ESC","EpiLC","E4.5_Epiblast","E5.5_Epiblast","E6.5_Epiblast","E6.5_Primitive_Streak",
              "E7.5_Ectoderm","E7.5_Endoderm","E7.5_Mesoderm",
              "E11_Brain.fore","E11_Brain.mid","E11_Brain.hind","E11_Liver","E11_Heart","E11_Adipose")

print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())

tissue.plot=c("somatic_adrenal_gland","somatic_gastrocnemius","somatic_heart","somatic_left_cerebral_cortex","somatic_hippocampus")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())


##Somatic plots##


dev.off()



##Ploting paper fig 2 figures
pdf("./Plots/Gene_expression/Mouse_Prim_nonprim_comparisons_paper_format.pdf")

enh.plot=c("mECT_prim")
sub.df=enh.coding[enh.coding$Enh.class%in%enh.plot,]
tissue.plot=c("E4.5_Epiblast",
              "E7.5_Ectoderm","E7.5_Endoderm","E7.5_Mesoderm")
sub.df=enh.coding[enh.coding$Tissue%in%tissue.plot,]
print(ggplot(sub.df[sub.df$Enh.class%in%c("mECT_prim"),],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Enh.class%in%c("mECT_nonprim"),],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Enh.class%in%c("mEND_prim"),],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Enh.class%in%c("mEND_nonprim"),],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Enh.class%in%c("mMES_prim"),],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Enh.class%in%c("mMES_nonprim"),],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())

tissue.plot=c("E11_Brain.fore","E11_Liver","E11_Heart","somatic_adrenal_gland","somatic_left_cerebral_cortex","somatic_heart")
sub.df=enh.coding[enh.coding$Tissue%in%tissue.plot,]
print(ggplot(sub.df[sub.df$Enh.class%in%c("mECT_prim","mECT_nonprim"),],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Enh.class%in%c("mEND_prim","mEND_nonprim"),],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df[sub.df$Enh.class%in%c("mMES_prim","mMES_nonprim"),],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())

dev.off()



#Aim to generate a list of Non-primed and Primed (outside of tissue) plots
nonprim.classes=gsub("_prim","_nonprim",prim.classes)

tissues=unique(enh.coding$Tissue)
tissue.prim.lookup=list(
  "ESC"=prim.classes[!prim.classes%in%c("Embryonic")],
  "EpiLC"=prim.classes[!prim.classes%in%c("Embryonic")],
  "E4.5_Epiblast"=prim.classes[!prim.classes%in%c("Embryonic")],
  "E5.5_Epiblast"=prim.classes[!prim.classes%in%c("Embryonic")],
  "E6.5_Epiblast"=prim.classes[!prim.classes%in%c("Embryonic")],
  "E6.5_Primitive_Streak"=prim.classes[!prim.classes%in%c("Embryonic")],
  "E7.5_Ectoderm"=prim.classes[!prim.classes%in%c("mECT_prim")],
  "E7.5_Endoderm"=prim.classes[!prim.classes%in%c("mEND_prim")],
  "E7.5_Mesoderm"=prim.classes[!prim.classes%in%c("mMES_prim")],
  "E11_Brain.fore"=prim.classes[!prim.classes%in%c("forebrain_prim","midbrain_prim","hindbrain_prim")],
  "E11_Brain.mid"=prim.classes[!prim.classes%in%c("forebrain_prim","midbrain_prim","hindbrain_prim")],
  "E11_Brain.hind"=prim.classes[!prim.classes%in%c("forebrain_prim","midbrain_prim","hindbrain_prim")],
  "E11_Liver"=prim.classes[!prim.classes%in%c("liver_prim")],
  "E11_Heart"=prim.classes[!prim.classes%in%c("heart_prim")],
  "E11_Adipose"=prim.classes[!prim.classes%in%c()],
  "somatic_adrenal_gland"=prim.classes[!prim.classes%in%c()],
  "somatic_gastrocnemius"=prim.classes[!prim.classes%in%c()],
  "somatic_heart"=prim.classes[!prim.classes%in%c()],
  "somatic_left_cerebral_cortex"=prim.classes[!prim.classes%in%c()],
  "somatic_hippocampus"=prim.classes[!prim.classes%in%c()])

tissue.nonprim.lookup=list(
  "ESC"=nonprim.classes[!nonprim.classes%in%c("Embryonic")],
  "EpiLC"=nonprim.classes[!nonprim.classes%in%c("Embryonic")],
  "E4.5_Epiblast"=nonprim.classes[!nonprim.classes%in%c("Embryonic")],
  "E5.5_Epiblast"=nonprim.classes[!nonprim.classes%in%c("Embryonic")],
  "E6.5_Epiblast"=nonprim.classes[!nonprim.classes%in%c("Embryonic")],
  "E6.5_Primitive_Streak"=nonprim.classes[!nonprim.classes%in%c("Embryonic")],
  "E7.5_Ectoderm"=nonprim.classes[!nonprim.classes%in%c("mECT_nonprim")],
  "E7.5_Endoderm"=nonprim.classes[!nonprim.classes%in%c("mEND_nonprim")],
  "E7.5_Mesoderm"=nonprim.classes[!nonprim.classes%in%c("mMES_nonprim")],
  "E11_Brain.fore"=nonprim.classes[!nonprim.classes%in%c("forebrain_nonprim","midbrain_nonprim","hindbrain_nonprim")],
  "E11_Brain.mid"=nonprim.classes[!nonprim.classes%in%c("forebrain_nonprim","midbrain_nonprim","hindbrain_nonprim")],
  "E11_Brain.hind"=nonprim.classes[!nonprim.classes%in%c("forebrain_nonprim","midbrain_nonprim","hindbrain_nonprim")],
  "E11_Liver"=nonprim.classes[!nonprim.classes%in%c("liver_nonprim")],
  "E11_Heart"=nonprim.classes[!nonprim.classes%in%c("heart_nonprim")],
  "E11_Adipose"=nonprim.classes[!nonprim.classes%in%c()],
  "somatic_adrenal_gland"=nonprim.classes[!nonprim.classes%in%c()],
  "somatic_gastrocnemius"=nonprim.classes[!nonprim.classes%in%c()],
  "somatic_heart"=nonprim.classes[!nonprim.classes%in%c()],
  "somatic_left_cerebral_cortex"=nonprim.classes[!nonprim.classes%in%c()],
  "somatic_hippocampus"=nonprim.classes[!nonprim.classes%in%c()])

prim.rows=list()
nonprim.rows=list()
for(i in tissues){
  df.sub=enh.coding[enh.coding$Tissue==i,]
  prim.sub=df.sub[df.sub$Enh.class%in%tissue.prim.lookup[[i]],]
  prim.sub$Enh.class=rep("Other_Primed")
  prim.sub=unique(prim.sub)
  non.sub=df.sub[df.sub$Enh.class%in%tissue.nonprim.lookup[[i]],]
  non.sub$Enh.class=rep("Other_Nonprimed")
  non.sub=unique(non.sub)
  prim.rows[[i]]=prim.sub
  nonprim.rows[[i]]=non.sub
  print(i)
}

tog.df=rbind(do.call(rbind,prim.rows),do.call(rbind,nonprim.rows))
tog.df$Tissue=factor(tog.df$Tissue,levels=plot.order)
tissue.plot=c("ESC","EpiLC","E4.5_Epiblast","E5.5_Epiblast","E6.5_Epiblast","E6.5_Primitive_Streak",
              "E7.5_Ectoderm","E7.5_Endoderm","E7.5_Mesoderm",
              "E11_Brain.fore","E11_Brain.mid","E11_Brain.hind","E11_Liver","E11_Heart","E11_Adipose",
              "somatic_adrenal_gland","somatic_gastrocnemius","somatic_heart","somatic_left_cerebral_cortex","somatic_hippocampus")
tog.df$Tissue=factor(tog.df$Tissue,levels = tissue.plot)

pdf("./Plots/Gene_expression/Mouse_Primed_vs_nonprimed_in_other_tissues.pdf")
p=ggboxplot(tog.df,x="Tissue",y="percentile",color="Enh.class")
p+stat_compare_means(aes(group=Enh.class),label="p.signif",method="wilcox.test")+stat_compare_means(aes(group=Enh.class),method="wilcox.test",label.y=120)
dev.off()



####Human####

setwd("G:/Headstone/Enh_priming/")

library(ggplot2)
library(tidyverse)
library(qlcMatrix)
library(ggpubr)


##Loading protein coding genes

tss=read.delim("./Annotations/Human/ENSEMBL_coding_gene_TSS.txt",h=F)
gene.ids=unique(tss[,4])


##Function to add coding genes which were not covered by a dataset and sets them to the minimal observed value
add.not.in=function(df){
  not.in.ids=gene.ids[!gene.ids%in%df[,1]]
  min.val=min(df[,2:ncol(df)]) 
  not.in.mat=as.data.frame(matrix(data=rep(min.val),ncol = ncol(df),nrow=length(not.in.ids)))
  not.in.mat[,1]=not.in.ids
  colnames(not.in.mat)=colnames(df)
  return(rbind(df,not.in.mat))
}


##Loading and formating all relevant expression datasets
xie=read.delim("./Data/human/RNA/Xie_et_al_RNA_logRPKM_all_types.txt")[,c(7,13:16)]
xie=add.not.in(xie)
xie=pivot_longer(xie,cols = -1, names_to = "Tissue",values_to = "Expr")

fetal=read.delim("./Data/human/RNA/Fetal_Atlas_GSE156793_S5_gene_fraction_tissue.txt",sep = ",")
fetal$RowID=unlist(lapply(fetal$RowID,function(x){unlist(strsplit(x,"\\."))[1]}))
colnames(fetal)[2:ncol(fetal)]=paste0("Fetal_",colnames(fetal)[2:ncol(fetal)])
colnames(fetal)[1]="ID"
fetal=add.not.in(fetal)
fetal=pivot_longer(fetal,cols = -1, names_to = "Tissue",values_to = "Expr")

somatic=read.delim("./Data/human/RNA/Cui_Somatic_RNA_logRPKM_all_types.txt")[,c(7,13:31)]
colnames(somatic)[2:20]=paste0("Somatic_",colnames(somatic)[2:20])
somatic=add.not.in(somatic)
somatic=pivot_longer(somatic,cols = -1, names_to = "Tissue",values_to = "Expr")
#combining all datasets
data=do.call(rbind,list(xie,fetal,somatic))
#subsetting to only protein coding genes
data.coding=data[data$ID%in%gene.ids,]
#calculating percentile positioning values
data.coding=data.coding%>%
  group_by(Tissue)%>% # group by lo and li
  mutate(percentile=percent_rank(Expr)*100) # make new column

fetal.order=c("Fetal_Cerebellum","Fetal_Cerebrum" ,"Fetal_Eye" ,
              "Fetal_Adrenal","Fetal_Muscle","Fetal_Kidney",
              "Fetal_Spleen" ,"Fetal_Thymus","Fetal_Heart",
              "Fetal_Intestine","Fetal_Liver","Fetal_Lung",
              "Fetal_Pancreas","Fetal_Stomach","Fetal_Placenta")
somatic.order = c("Somatic_brain",  "Somatic_hypothalamus",  "Somatic_skin",
                  "Somatic_breast","Somatic_bonemarrow","Somatic_kidney",
                  "Somatic_lymphnodes","Somatic_ovary","Somatic_testis",
                  "Somatic_bladder","Somatic_stomach","Somatic_colonsigmoid","Somatic_colontransverse",
                  "Somatic_liver","Somatic_heart","Somatic_lung" ,
                  "Somatic_pancreas" ,"Somatic_prostate" ,"Somatic_placenta")


tissue.plot.order=c("H1","NPC","ME","MSC",fetal.order,somatic.order)
tissue.sub=c("H1","NPC","ME","MSC","Fetal_Cerebellum","Fetal_Heart","Fetal_Liver","Somatic_brain","Somatic_heart","Somatic_liver")
fetal.sub=c("Fetal_Cerebellum","Fetal_Adrenal","Fetal_Heart","Fetal_Liver","Fetal_Lung","Fetal_Stomach","Fetal_Kidney")
data.coding$Tissue=factor(data.coding$Tissue,levels=tissue.plot.order)
colnames(data.coding)[1]="Gene.ID"
##Checking the normalisation
ggplot(data.coding,aes(x=Tissue,y=percentile))+geom_boxplot(outlier.shape = NA)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####Getting the human enhancer groups and gene associations####
enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(enh.lookup)=c("Enh.ID","Enh.class")
classes=unique(enh.lookup$Enh.class)

prox=read.delim("./Annotations/Human/Prox_genes/Nearest_expressed_gene_20kb.txt",h=F)
colnames(prox)=c("Enh.ID","Gene.ID")
hic=read.delim("./Annotations/Human/HiC_genes/hESC_HiC_genes.txt",h=F)
colnames(hic)=c("Enh.ID","Gene.ID")
comb=rbind(prox,hic)

comb=merge(comb,enh.lookup,by="Enh.ID")



embryonic.genes=unique(comb$Gene.ID[grepl("hESC_specific",comb$Enh.class)])


rel.classes=c("NPC_prim"  , "NPC_nonprim","ME_prim" ,"ME_nonprim","MSC_prim", "MSC_nonprim","NEC_prim","NEC_nonprim",
              "7W_brain_prim"  , "7W_brain_nonprim" ,"Fetal_adrenal_prim","Fetal_adrenal_nonprim","Fetal_kidney_prim" ,
              "Fetal_kidney_nonprim","Fetal_liver_prim","Fetal_liver_nonprim", "Fetal_lung_prim"   , "Fetal_lung_nonprim",
              "Fetal_stomach_specific", "Fetal_stomach_prim", "Fetal_stomach_nonprim","Fetal_ventricle_prim"  ,"Fetal_ventricle_nonprim")

prim.classes=rel.classes[grepl("_prim",rel.classes)]
df.rows=list()
for(i in rel.classes){
  if(grepl("_prim",i)){
    gene.list=unique(comb$Gene.ID[comb$Enh.class==i])
    other.prim.classes=prim.classes[prim.classes!=i]
    other.prim=comb$Gene.ID[comb$Enh.class%in%other.prim.classes]
    df.rows[[paste0(i,".embryonic")]]=data.frame("Enh.class"=rep(paste0(i,".embryonic")),"Gene.ID"=gene.list[gene.list%in%embryonic.genes&!gene.list%in%other.prim])
    gene.list=gene.list[!gene.list%in%embryonic.genes]
    other.prim.classes=prim.classes[prim.classes!=i]
    other.prim=comb$Gene.ID[comb$Enh.class%in%other.prim.classes]
    df.rows[[paste0(i,".specific")]]=data.frame("Enh.class"=rep(i),"Gene.ID"=gene.list[!gene.list%in%other.prim])
  }else   if(grepl("_nonprim",i)){
    gene.list=unique(comb$Gene.ID[comb$Enh.class==i])
    #gene.list=gene.list[!gene.list%in%embryonic.genes]
    prim.list=comb$Gene.ID[grepl("_prim",comb$Enh.class)]
    gene.list=gene.list[!gene.list%in%prim.list]
    df.rows[[i]]=data.frame("Enh.class"=rep(i),"Gene.ID"=gene.list)
  }
}

df.tog=do.call(rbind,df.rows)

enh.comb=merge(data.coding,df.tog,by="Gene.ID")

embryonic.df=data.coding[data.coding$Gene.ID%in%gene.ids,]
embryonic.df$Enh.class=rep("hESC")
data.comb=rbind(enh.comb,embryonic.df[,colnames(enh.comb)])
data.comb$Tissue=factor(data.comb$Tissue,levels=tissue.plot.order)


##Ploting final comparisons
pdf("./Plots/Gene_expression/Human_Prim_nonprim_comparisons_across_dev.pdf")
enh.plot=c("NPC_prim","NPC_nonprim","ME_prim","ME_nonprim","MSC_prim","MSC_nonprim")
sub.df=data.comb[data.comb$Enh.class%in%enh.plot,]
tissue.plot=c("H1","NPC","ME","MSC")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
tissue.plot=c("H1","NPC","ME","MSC","Fetal_Liver","Fetal_Lung","Fetal_Cerebellum","Somatic_liver","Somatic_lung","Somatic_brain")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())

enh.plot=c("7W_brain_prim","7W_brain_nonprim")
sub.df=data.comb[data.comb$Enh.class%in%enh.plot,]
sub.df$Enh.class=factor(sub.df$Enh.class,levels=enh.plot)
tissue.plot=c("H1","NPC","ME","MSC","Fetal_Liver","Fetal_Lung","Fetal_Cerebellum")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
tissue.plot=c("H1","NPC","ME","MSC","Fetal_Liver","Fetal_Lung","Fetal_Cerebellum","Somatic_liver","Somatic_lung","Somatic_brain")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())

enh.plot=c("Fetal_adrenal_prim","Fetal_adrenal_nonprim","Fetal_kidney_prim","Fetal_kidney_nonprim","Fetal_liver_prim","Fetal_liver_nonprim")
sub.df=data.comb[data.comb$Enh.class%in%enh.plot,]
sub.df$Enh.class=factor(sub.df$Enh.class,levels=enh.plot)
tissue.plot=c("H1","NPC","ME","MSC","Fetal_Adrenal","Fetal_Kidney","Fetal_Liver")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
tissue.plot=c("H1","NPC","ME","MSC","Fetal_Adrenal","Fetal_Kidney","Fetal_Liver","Fetal_Liver","Fetal_Lung","Fetal_Cerebellum","Somatic_liver","Somatic_lung","Somatic_brain")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())

enh.plot=c("Fetal_lung_prim","Fetal_lung_nonprim","Fetal_stomach_prim","Fetal_stomach_nonprim","Fetal_ventricle_prim","Fetal_ventricle_nonprim")
sub.df=data.comb[data.comb$Enh.class%in%enh.plot,]
sub.df$Enh.class=factor(sub.df$Enh.class,levels=enh.plot)
tissue.plot=c("H1","NPC","ME","MSC","Fetal_Lung","Fetal_Stomach","Fetal_Heart")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
tissue.plot=c("H1","NPC","ME","MSC","Fetal_Lung","Fetal_Stomach","Fetal_Heart","Fetal_Liver","Fetal_Lung","Fetal_Cerebellum","Somatic_liver","Somatic_lung","Somatic_brain")
print(ggplot(sub.df[sub.df$Tissue%in%tissue.plot,],aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())

dev.off()


pdf("./Plots/Gene_expression/Human_Prim_nonprim_ALL_comparisons.pdf")
enh.plot=c("NPC_prim","NPC_nonprim","ME_prim","ME_nonprim","MSC_prim","MSC_nonprim")
sub.df=data.comb[data.comb$Enh.class%in%enh.plot,]
print(ggplot(sub.df,aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df,aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())

enh.plot=c("7W_brain_prim","7W_brain_nonprim")
sub.df=data.comb[data.comb$Enh.class%in%enh.plot,]
print(ggplot(sub.df,aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df,aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())

enh.plot=c("Fetal_adrenal_prim","Fetal_adrenal_nonprim","Fetal_kidney_prim","Fetal_kidney_nonprim","Fetal_liver_prim","Fetal_liver_nonprim")
sub.df=data.comb[data.comb$Enh.class%in%enh.plot,]
print(ggplot(sub.df,aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df,aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())
enh.plot=c("Fetal_lung_prim","Fetal_lung_nonprim","Fetal_stomach_prim","Fetal_stomach_nonprim","Fetal_ventricle_prim","Fetal_ventricle_nonprim")
sub.df=data.comb[data.comb$Enh.class%in%enh.plot,]
print(ggplot(sub.df,aes(x=Tissue,y=percentile,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+theme_bw())
print(ggplot(sub.df,aes(x=Enh.class,y=percentile,fill=Tissue))+geom_boxplot(outlier.shape = NA)+theme_bw())

dev.off()






prim.classes=rel.classes[grepl("_prim",rel.classes)]
prim.classes=prim.classes[grepl("NPC|ME|MSC",prim.classes)]
#prim.classes=prim.classes[grepl("Fetal",prim.classes)]

nonprim.classes=gsub("_prim","_nonprim",prim.classes)

tissues=unique(gast.comb$Tissue)
tissue.prim.lookup=list(
  "H1"=prim.classes[!prim.classes%in%c("hESC_specific")],
  "NPC"=prim.classes[!prim.classes%in%c("NPC_prim")],
  "ME"=prim.classes[!prim.classes%in%c("ME_prim")],
  "MSC"=prim.classes[!prim.classes%in%c("MSC_prim")],
  "Fetal_Cerebellum"=prim.classes[!prim.classes%in%c("7W_brain_prim")],
  "Fetal_Cerebrum"=prim.classes[!prim.classes%in%c("7W_brain_prim")],
  "Fetal_Eye"=prim.classes[!prim.classes%in%c("")],
  "Fetal_Adrenal"=prim.classes[!prim.classes%in%c("Fetal_adrenal_prim")],
  "Fetal_Muscle"=prim.classes[!prim.classes%in%c("")],
  "Fetal_Kidney"=prim.classes[!prim.classes%in%c("Fetal_kidney_prim")],
  "Fetal_Spleen"=prim.classes[!prim.classes%in%c("")],
  "Fetal_Thymus"=prim.classes[!prim.classes%in%c("")],
  "Fetal_Heart"=prim.classes[!prim.classes%in%c("")],
  "Fetal_Intestine"=prim.classes[!prim.classes%in%c("")],
  "Fetal_Liver"=prim.classes[!prim.classes%in%c("Fetal_liver_prim")],
  "Fetal_Lung"=prim.classes[!prim.classes%in%c("Fetal_lung_prim")],
  "Fetal_Pancreas"=prim.classes[!prim.classes%in%c()],
  "Fetal_Stomach"=prim.classes[!prim.classes%in%c("Fetal_stomach_prim")],
  "Fetal_Placenta"=prim.classes[!prim.classes%in%c()],
  "Somatic_brain"=prim.classes[!prim.classes%in%c()],
  "Somatic_hypothalamus"=prim.classes[!prim.classes%in%c()],
  "Somatic_skin"=prim.classes[!prim.classes%in%c()],
  "Somatic_breast"=prim.classes[!prim.classes%in%c()],
  "Somatic_bonemarrow"=prim.classes[!prim.classes%in%c()],
  "Somatic_kidney"=prim.classes[!prim.classes%in%c()],
  "Somatic_lymphnodes"=prim.classes[!prim.classes%in%c()],
  "Somatic_ovary"=prim.classes[!prim.classes%in%c()],
  "Somatic_testis"=prim.classes[!prim.classes%in%c()],
  "Somatic_bladder"=prim.classes[!prim.classes%in%c()],
  "Somatic_stomach"=prim.classes[!prim.classes%in%c()],
  "Somatic_colonsigmoid"=prim.classes[!prim.classes%in%c()],
  "Somatic_colontransverse"=prim.classes[!prim.classes%in%c()],
  "Somatic_liver"=prim.classes[!prim.classes%in%c()],
  "Somatic_heart"=prim.classes[!prim.classes%in%c()],
  "Somatic_lung"=prim.classes[!prim.classes%in%c()],
  "Somatic_pancreas"=prim.classes[!prim.classes%in%c()],
  "Somatic_prostate"=prim.classes[!prim.classes%in%c()],
  "Somatic_placenta"=prim.classes[!prim.classes%in%c()])
tissue.nonprim.lookup=list(
  "H1"=nonprim.classes[!nonprim.classes%in%c("hESC_specific")],
  "NPC"=nonprim.classes[!nonprim.classes%in%c("NPC_nonprim")],
  "ME"=nonprim.classes[!nonprim.classes%in%c("ME_nonprim")],
  "MSC"=nonprim.classes[!nonprim.classes%in%c("MSC_nonprim")],
  "Fetal_Cerebellum"=nonprim.classes[!nonprim.classes%in%c("7W_brain_nonprim")],
  "Fetal_Cerebrum"=nonprim.classes[!nonprim.classes%in%c("7W_brain_nonprim")],
  "Fetal_Eye"=nonprim.classes[!nonprim.classes%in%c("")],
  "Fetal_Adrenal"=nonprim.classes[!nonprim.classes%in%c("Fetal_adrenal_nonprim")],
  "Fetal_Muscle"=nonprim.classes[!nonprim.classes%in%c("")],
  "Fetal_Kidney"=nonprim.classes[!nonprim.classes%in%c("Fetal_kidney_nonprim")],
  "Fetal_Spleen"=nonprim.classes[!nonprim.classes%in%c("")],
  "Fetal_Thymus"=nonprim.classes[!nonprim.classes%in%c("")],
  "Fetal_Heart"=nonprim.classes[!nonprim.classes%in%c("")],
  "Fetal_Intestine"=nonprim.classes[!nonprim.classes%in%c("")],
  "Fetal_Liver"=nonprim.classes[!nonprim.classes%in%c("Fetal_liver_nonprim")],
  "Fetal_Lung"=nonprim.classes[!nonprim.classes%in%c("Fetal_lung_nonprim")],
  "Fetal_Pancreas"=nonprim.classes[!nonprim.classes%in%c()],
  "Fetal_Stomach"=nonprim.classes[!nonprim.classes%in%c("Fetal_stomach_nonprim")],
  "Fetal_Placenta"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_brain"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_hypothalamus"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_skin"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_breast"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_bonemarrow"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_kidney"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_lymphnodes"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_ovary"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_testis"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_bladder"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_stomach"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_colonsigmoid"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_colontransverse"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_liver"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_heart"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_lung"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_pancreas"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_prostate"=nonprim.classes[!nonprim.classes%in%c()],
  "Somatic_placenta"=nonprim.classes[!nonprim.classes%in%c()])


prim.rows=list()
nonprim.rows=list()
for(i in tissues){
  df.sub=data.comb[data.comb$Tissue==i,]
  prim.sub=df.sub[df.sub$Enh.class%in%tissue.prim.lookup[[i]],]
  prim.sub$Enh.class=rep("Other_Primed")
  prim.sub=unique(prim.sub)
  non.sub=df.sub[df.sub$Enh.class%in%tissue.nonprim.lookup[[i]],]
  non.sub$Enh.class=rep("Other_Nonprimed")
  non.sub=unique(non.sub)
  prim.rows[[i]]=prim.sub
  nonprim.rows[[i]]=non.sub
  print(i)
}

tog.df=rbind(do.call(rbind,prim.rows),do.call(rbind,nonprim.rows))
tog.df$Tissue=factor(tog.df$Tissue,levels=tissue.plot.order)
tissue.plot=c("H1","NPC","ME","MSC","Fetal_Liver","Fetal_Lung","Fetal_Cerebellum","Somatic_liver","Somatic_lung","Somatic_brain")

pdf("./Plotsc/Human_other_nonprim_expr.pdf")
p=ggboxplot(tog.df[tog.df$Tissue%in%tissue.plot,],x="Tissue",y="percentile",color="Enh.class")
p+stat_compare_means(aes(group=Enh.class),label="p.signif",method="wilcox.test")+stat_compare_means(aes(group=Enh.class),method="wilcox.test",label.y=120)
p=ggboxplot(tog.df,x="Tissue",y="percentile",color="Enh.class")
p+stat_compare_means(aes(group=Enh.class),label="p.signif",method="wilcox.test")+stat_compare_means(aes(group=Enh.class),method="wilcox.test",label.y=120)
dev.off()
