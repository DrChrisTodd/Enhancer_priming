setwd("G:/Headstone/Enh_priming/")

library(eulerr)

####Generating a TF ortho lookup file####
tfs=read.delim("./Annotations/Human/TFs.txt",h=F)[,1]
ortho=read.delim("./Rebuttal/Annotations/Human_mouse_orthologous_biomart.txt")[,c(6,1,4,3,5)]
h.gene.names=read.delim("./Annotations/Human/ENSEMBL_IDs_to_gene_names.txt")
m.gene.names=read.delim("./Annotations/Mouse/ENSEMBL_IDs_to_gene_names.txt")

h.orth=ortho[ortho$Gene.name%in%tfs,]
m.ortho=ortho[toupper(ortho$Mouse.gene.name)%in%tfs,]
tf.ortho=rbind(h.orth,m.ortho)
tf.ortho=unique(tf.ortho)

df=data.frame(tfs,"in.hum"=rep(NA),"in.mouse"=rep(NA),"in.h.ortho"=rep(NA),"in.m.ortho"=rep(NA))
df$in.hum=df$tfs%in%h.gene.names$Gene.name
df$in.mouse=df$tfs%in%toupper(m.gene.names$Gene.name)
df$in.h.ortho=df$tfs%in%ortho$Gene.name
df$in.m.ortho=df$tfs%in%toupper(ortho$Mouse.gene.name)

h.only.tf=df[df$in.hum&!df$in.mouse,]
m.only.tf=df[df$in.mouse&!df$in.hum,]
h.only.ortho=h.only.tf[h.only.tf$in.h.ortho,]
h.tru.only=h.only.tf[!h.only.tf$in.h.ortho,]
m.tru.only=m.only.tf[!m.only.tf$in.m.ortho,]

h.sub=h.gene.names[h.gene.names$Gene.name%in%h.tru.only$tfs,]
h.df=data.frame("Gene.name"=h.sub$Gene.name,"Gene.stable.ID"=h.sub$Gene.stable.ID,
                "Mouse.gene.name"=rep(NA),"Mouse.gene.stable.ID"=rep(NA),
                "Mouse.homology.type"=rep("none"))
m.sub=unique(m.gene.names[toupper(m.gene.names$Gene.name)%in%m.tru.only$tfs,1:2])
m.df=data.frame("Gene.name"=rep(NA),"Gene.stable.ID"=rep(NA),
                "Mouse.gene.name"=m.sub$Gene.name,"Mouse.gene.stable.ID"=m.sub$Gene.stable.ID,
                "Mouse.homology.type"=rep("none"))
tf.tog=do.call(rbind,list(tf.ortho,h.df,m.df))
tf.tog$TF.ID=paste0("TF_",seq(1:nrow(tf.tog)))
write.table(tf.tog,"./Rebuttal/Annotations/TF_ortho_lookup.txt",sep="\t",col.names = T,row.names = F,quote = F)

####Checking for shared TFs####
tfs=read.delim("./Rebuttal/Annotations/TF_ortho_lookup.txt")
h.gene.names=read.delim("./Annotations/Human/ENSEMBL_IDs_to_gene_names.txt")
m.gene.names=read.delim("./Annotations/Mouse/ENSEMBL_IDs_to_gene_names.txt")
ortho=read.delim("./Rebuttal/Annotations/Human_mouse_orthologous_biomart.txt")[,c(6,1,4,3,5)]
ortho$Ortho.ID=paste0("Ortho_",seq(1:nrow(ortho)))

{####Loading Mouse Enhancer-gene associations####
  
  enh.lookup=read.delim("./Annotations/Mouse/Comb_Enh_class_lookup.txt",h=F)
  colnames(enh.lookup)=c("ID","Enh.class")
  classes=unique(enh.lookup$Enh.class)
  
  epi.hic=read.delim("./Annotations/Mouse/HiC_genes/EpiLC_HiC_genes.txt",h=F)
  prox=read.delim("./Annotations/Mouse/Prox_genes/Nearest_expressed_gene_20kb.txt",h=F)
  colnames(epi.hic)=colnames(prox)=c("ID","Gene.ID")
  
  m.comb=rbind(epi.hic,prox)
  m.comb=merge(m.comb,enh.lookup,by="ID")
  
  ###Treating the enhancer-gene associations
  ##Removing Embryonic enhancer associated genes, and genes regulated by more than one lineage of enhancers
  embryonic.genes=unique(m.comb$Gene.ID[grepl("ESC|EpiLC",m.comb$Enh.class)])
  
  rel.classes=c("mECT_prim","mECT_nonprim","mEND_prim","mEND_nonprim","mMES_prim","mMES_nonprim","liver_prim","liver_nonprim",    
                "heart_prim","heart_nonprim","forebrain_prim","forebrain_nonprim","midbrain_prim","midbrain_nonprim","hindbrain_prim","hindbrain_nonprim")
  prim.classes=rel.classes[grepl("_prim",rel.classes)]
  df.rows=list()
  for(i in rel.classes){
    if(grepl("_prim",i)){
      gene.list=unique(m.comb$Gene.ID[m.comb$Enh.class==i])
      other.prim.classes=prim.classes[prim.classes!=i]
      other.prim=m.comb$Gene.ID[m.comb$Enh.class%in%other.prim.classes]
      df.rows[[paste0(i,".embryonic")]]=data.frame("Enh.class"=rep(paste0(i,".embryonic")),"Gene.ID"=gene.list[gene.list%in%embryonic.genes&!gene.list%in%other.prim])
      gene.list=gene.list[!gene.list%in%embryonic.genes]
      other.prim.classes=prim.classes[prim.classes!=i]
      other.prim=m.comb$Gene.ID[m.comb$Enh.class%in%other.prim.classes]
      df.rows[[paste0(i,".specific")]]=data.frame("Enh.class"=rep(i),"Gene.ID"=gene.list[!gene.list%in%other.prim])
    }else   if(grepl("_nonprim",i)){
      gene.list=unique(m.comb$Gene.ID[m.comb$Enh.class==i])
      # gene.list=gene.list[!gene.list%in%embryonic.genes]
      prim.list=m.comb$Gene.ID[grepl("_prim",m.comb$Enh.class)]
      gene.list=gene.list[!gene.list%in%prim.list]
      df.rows[[i]]=data.frame("Enh.class"=rep(i),"Gene.ID"=gene.list)
    }
  }
  df.tog=do.call(rbind,df.rows)
}
{####Getting the human enhancer groups and gene associations####
  enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
  colnames(enh.lookup)=c("Enh.ID","Enh.class")
  classes=unique(enh.lookup$Enh.class)
  
  prox=read.delim("./Annotations/Human/Prox_genes/Nearest_expressed_gene_20kb.txt",h=F)
  colnames(prox)=c("Enh.ID","Gene.ID")
  hic=read.delim("./Annotations/Human/HiC_genes/hESC_HiC_genes.txt",h=F)
  colnames(hic)=c("Enh.ID","Gene.ID")
  h.comb=rbind(prox,hic)
  
  h.comb=merge(h.comb,enh.lookup,by="Enh.ID")
  
}

head(h.comb)

coi="hNPC_vs_mECT"
h.class="NPC_prim"
m.class="mECT_prim"

{
h.sub=h.comb[h.comb$Enh.class==h.class,]
m.sub=m.comb[m.comb$Enh.class==m.class,]

h.sub.ortho=h.sub[h.sub$Gene.ID%in%ortho$Gene.stable.ID,]
m.sub.ortho=m.sub[m.sub$Gene.ID%in%ortho$Mouse.gene.stable.ID,]

#Getting venn data
all.human=unique(h.sub$Gene.ID)
h.not.ortho=all.human[!all.human%in%ortho$Gene.stable.ID]
human.ortho=ortho$Ortho.ID[ortho$Gene.stable.ID%in%all.human]
mouse.all.ids=unique(m.sub$Gene.ID)
m.not.ortho=mouse.all.ids[!mouse.all.ids%in%ortho$Mouse.gene.stable.ID]
mouse.ortho=ortho$Ortho.ID[ortho$Mouse.gene.stable.ID%in%mouse.all.ids]
all.gene.ids=unique(c(h.not.ortho,human.ortho,m.not.ortho,mouse.ortho))

gene.mat=matrix(nrow=length(all.gene.ids),ncol=4)
colnames(gene.mat)=c("All_Human","Ortho_Human","Ortho_Mouse","All_Mouse")
gene.mat[,1]=all.gene.ids%in%c(h.not.ortho,human.ortho)
gene.mat[,2]=all.gene.ids%in%human.ortho
gene.mat[,3]=all.gene.ids%in%mouse.ortho
gene.mat[,4]=all.gene.ids%in%c(m.not.ortho,mouse.ortho)
pdf(paste0("./Rebuttal/Plots/Conserved_genes/",coi,"_num_shared_genes.pdf"))
fit=euler(gene.mat)
plot(fit,quantities = T)
dev.off()
shared.genes=ortho$Gene.name[ortho$Ortho.ID%in%human.ortho[human.ortho%in%mouse.ortho]]
write.table(shared.genes,paste0("./Rebuttal/Plots/Conserved_Genes/",coi,"_shared_genes_ids.txt"),col.names = F,row.names = F,quote = F,sep = "\t")
}

coi="hBrain_vs_mBrain"
h.class="7W_brain_prim"
m.class="brain_prim"

{
  h.sub=h.comb[h.comb$Enh.class==h.class,]
  m.sub=m.comb[m.comb$Enh.class==m.class,]
  
  h.sub.ortho=h.sub[h.sub$Gene.ID%in%ortho$Gene.stable.ID,]
  m.sub.ortho=m.sub[m.sub$Gene.ID%in%ortho$Mouse.gene.stable.ID,]
  
  #Getting venn data
  all.human=unique(h.sub$Gene.ID)
  h.not.ortho=all.human[!all.human%in%ortho$Gene.stable.ID]
  human.ortho=ortho$Ortho.ID[ortho$Gene.stable.ID%in%all.human]
  mouse.all.ids=unique(m.sub$Gene.ID)
  m.not.ortho=mouse.all.ids[!mouse.all.ids%in%ortho$Mouse.gene.stable.ID]
  mouse.ortho=ortho$Ortho.ID[ortho$Mouse.gene.stable.ID%in%mouse.all.ids]
  all.gene.ids=unique(c(h.not.ortho,human.ortho,m.not.ortho,mouse.ortho))
  
  gene.mat=matrix(nrow=length(all.gene.ids),ncol=4)
  colnames(gene.mat)=c("All_Human","Ortho_Human","Ortho_Mouse","All_Mouse")
  gene.mat[,1]=all.gene.ids%in%c(h.not.ortho,human.ortho)
  gene.mat[,2]=all.gene.ids%in%human.ortho
  gene.mat[,3]=all.gene.ids%in%mouse.ortho
  gene.mat[,4]=all.gene.ids%in%c(m.not.ortho,mouse.ortho)
  pdf(paste0("./Rebuttal/Plots/Conserved_genes/",coi,"_num_shared_genes.pdf"))
  fit=euler(gene.mat)
  plot(fit,quantities = T)
  dev.off()
  shared.genes=ortho$Gene.name[ortho$Ortho.ID%in%human.ortho[human.ortho%in%mouse.ortho]]
  write.table(shared.genes,paste0("./Rebuttal/Plots/Conserved_Genes/",coi,"_shared_genes_ids.txt"),col.names = F,row.names = F,quote = F,sep = "\t")
}

coi="hliver_vs_mliver"
h.class="Fetal_liver_prim"
m.class="liver_prim"

{
  h.sub=h.comb[h.comb$Enh.class==h.class,]
  m.sub=m.comb[m.comb$Enh.class==m.class,]
  
  h.sub.ortho=h.sub[h.sub$Gene.ID%in%ortho$Gene.stable.ID,]
  m.sub.ortho=m.sub[m.sub$Gene.ID%in%ortho$Mouse.gene.stable.ID,]
  
  #Getting venn data
  all.human=unique(h.sub$Gene.ID)
  h.not.ortho=all.human[!all.human%in%ortho$Gene.stable.ID]
  human.ortho=ortho$Ortho.ID[ortho$Gene.stable.ID%in%all.human]
  mouse.all.ids=unique(m.sub$Gene.ID)
  m.not.ortho=mouse.all.ids[!mouse.all.ids%in%ortho$Mouse.gene.stable.ID]
  mouse.ortho=ortho$Ortho.ID[ortho$Mouse.gene.stable.ID%in%mouse.all.ids]
  all.gene.ids=unique(c(h.not.ortho,human.ortho,m.not.ortho,mouse.ortho))
  
  gene.mat=matrix(nrow=length(all.gene.ids),ncol=4)
  colnames(gene.mat)=c("All_Human","Ortho_Human","Ortho_Mouse","All_Mouse")
  gene.mat[,1]=all.gene.ids%in%c(h.not.ortho,human.ortho)
  gene.mat[,2]=all.gene.ids%in%human.ortho
  gene.mat[,3]=all.gene.ids%in%mouse.ortho
  gene.mat[,4]=all.gene.ids%in%c(m.not.ortho,mouse.ortho)
  pdf(paste0("./Rebuttal/Plots/Conserved_genes/",coi,"_num_shared_genes.pdf"))
  fit=euler(gene.mat)
  plot(fit,quantities = T)
  dev.off()
  shared.genes=ortho$Gene.name[ortho$Ortho.ID%in%human.ortho[human.ortho%in%mouse.ortho]]
  write.table(shared.genes,paste0("./Rebuttal/Plots/Conserved_Genes/",coi,"_shared_genes_ids.txt"),col.names = F,row.names = F,quote = F,sep = "\t")
}

coi="hNPC_vs_mECT_TFs"
h.class="NPC_prim"
m.class="mECT_prim"

{
  h.sub=h.comb[h.comb$Enh.class==h.class&h.comb$Gene.ID%in%tfs$Gene.stable.ID,]
  m.sub=m.comb[m.comb$Enh.class==m.class&m.comb$Gene.ID%in%tfs$Mouse.gene.stable.ID,]
  
  h.sub.ortho=h.sub[h.sub$Gene.ID%in%ortho$Gene.stable.ID,]
  m.sub.ortho=m.sub[m.sub$Gene.ID%in%ortho$Mouse.gene.stable.ID,]
  
  #Getting venn data
  all.human=unique(h.sub$Gene.ID)
  h.not.ortho=all.human[!all.human%in%ortho$Gene.stable.ID]
  human.ortho=ortho$Ortho.ID[ortho$Gene.stable.ID%in%all.human]
  mouse.all.ids=unique(m.sub$Gene.ID)
  m.not.ortho=mouse.all.ids[!mouse.all.ids%in%ortho$Mouse.gene.stable.ID]
  mouse.ortho=ortho$Ortho.ID[ortho$Mouse.gene.stable.ID%in%mouse.all.ids]
  all.gene.ids=unique(c(h.not.ortho,human.ortho,m.not.ortho,mouse.ortho))
  
  gene.mat=matrix(nrow=length(all.gene.ids),ncol=4)
  colnames(gene.mat)=c("All_Human","Ortho_Human","Ortho_Mouse","All_Mouse")
  gene.mat[,1]=all.gene.ids%in%c(h.not.ortho,human.ortho)
  gene.mat[,2]=all.gene.ids%in%human.ortho
  gene.mat[,3]=all.gene.ids%in%mouse.ortho
  gene.mat[,4]=all.gene.ids%in%c(m.not.ortho,mouse.ortho)
  pdf(paste0("./Rebuttal/Plots/Conserved_genes/",coi,"_num_shared_genes.pdf"))
  fit=euler(gene.mat)
  plot(fit,quantities = T)
  dev.off()
  shared.genes=ortho$Gene.name[ortho$Ortho.ID%in%human.ortho[human.ortho%in%mouse.ortho]]
  write.table(shared.genes,paste0("./Rebuttal/Plots/Conserved_Genes/",coi,"_shared_genes_ids.txt"),col.names = F,row.names = F,quote = F,sep = "\t")
}

coi="hBrain_vs_mBrain_TFs"
h.class="7W_brain_prim"
m.class="brain_prim"
{
  h.sub=h.comb[h.comb$Enh.class==h.class&h.comb$Gene.ID%in%tfs$Gene.stable.ID,]
  m.sub=m.comb[m.comb$Enh.class==m.class&m.comb$Gene.ID%in%tfs$Mouse.gene.stable.ID,]
  
  h.sub.ortho=h.sub[h.sub$Gene.ID%in%ortho$Gene.stable.ID,]
  m.sub.ortho=m.sub[m.sub$Gene.ID%in%ortho$Mouse.gene.stable.ID,]
  
  #Getting venn data
  all.human=unique(h.sub$Gene.ID)
  h.not.ortho=all.human[!all.human%in%ortho$Gene.stable.ID]
  human.ortho=ortho$Ortho.ID[ortho$Gene.stable.ID%in%all.human]
  mouse.all.ids=unique(m.sub$Gene.ID)
  m.not.ortho=mouse.all.ids[!mouse.all.ids%in%ortho$Mouse.gene.stable.ID]
  mouse.ortho=ortho$Ortho.ID[ortho$Mouse.gene.stable.ID%in%mouse.all.ids]
  all.gene.ids=unique(c(h.not.ortho,human.ortho,m.not.ortho,mouse.ortho))
  
  gene.mat=matrix(nrow=length(all.gene.ids),ncol=4)
  colnames(gene.mat)=c("All_Human","Ortho_Human","Ortho_Mouse","All_Mouse")
  gene.mat[,1]=all.gene.ids%in%c(h.not.ortho,human.ortho)
  gene.mat[,2]=all.gene.ids%in%human.ortho
  gene.mat[,3]=all.gene.ids%in%mouse.ortho
  gene.mat[,4]=all.gene.ids%in%c(m.not.ortho,mouse.ortho)
  pdf(paste0("./Rebuttal/Plots/Conserved_genes/",coi,"_num_shared_genes.pdf"))
  fit=euler(gene.mat)
  plot(fit,quantities = T)
  dev.off()
  shared.genes=ortho$Gene.name[ortho$Ortho.ID%in%human.ortho[human.ortho%in%mouse.ortho]]
  write.table(shared.genes,paste0("./Rebuttal/Plots/Conserved_Genes/",coi,"_shared_genes_ids.txt"),col.names = F,row.names = F,quote = F,sep = "\t")
}

coi="hliver_vs_mliver_TFs"
h.class="Fetal_liver_prim"
m.class="liver_prim"
{
  h.sub=h.comb[h.comb$Enh.class==h.class&h.comb$Gene.ID%in%tfs$Gene.stable.ID,]
  m.sub=m.comb[m.comb$Enh.class==m.class&m.comb$Gene.ID%in%tfs$Mouse.gene.stable.ID,]
  
  h.sub.ortho=h.sub[h.sub$Gene.ID%in%ortho$Gene.stable.ID,]
  m.sub.ortho=m.sub[m.sub$Gene.ID%in%ortho$Mouse.gene.stable.ID,]
  
  #Getting venn data
  all.human=unique(h.sub$Gene.ID)
  h.not.ortho=all.human[!all.human%in%ortho$Gene.stable.ID]
  human.ortho=ortho$Ortho.ID[ortho$Gene.stable.ID%in%all.human]
  mouse.all.ids=unique(m.sub$Gene.ID)
  m.not.ortho=mouse.all.ids[!mouse.all.ids%in%ortho$Mouse.gene.stable.ID]
  mouse.ortho=ortho$Ortho.ID[ortho$Mouse.gene.stable.ID%in%mouse.all.ids]
  all.gene.ids=unique(c(h.not.ortho,human.ortho,m.not.ortho,mouse.ortho))
  
  gene.mat=matrix(nrow=length(all.gene.ids),ncol=4)
  colnames(gene.mat)=c("All_Human","Ortho_Human","Ortho_Mouse","All_Mouse")
  gene.mat[,1]=all.gene.ids%in%c(h.not.ortho,human.ortho)
  gene.mat[,2]=all.gene.ids%in%human.ortho
  gene.mat[,3]=all.gene.ids%in%mouse.ortho
  gene.mat[,4]=all.gene.ids%in%c(m.not.ortho,mouse.ortho)
  pdf(paste0("./Rebuttal/Plots/Conserved_genes/",coi,"_num_shared_genes.pdf"))
  fit=euler(gene.mat)
  plot(fit,quantities = T)
  dev.off()
  shared.genes=ortho$Gene.name[ortho$Ortho.ID%in%human.ortho[human.ortho%in%mouse.ortho]]
  write.table(shared.genes,paste0("./Rebuttal/Plots/Conserved_Genes/",coi,"_shared_genes_ids.txt"),col.names = F,row.names = F,quote = F,sep = "\t")
}



summary(as.factor(h.TFs.comb$Enh.class))
writeLines(human.lookup$Gene.name[human.lookup$Gene.stable.ID%in%h.TFs.comb$Gene.ID[h.TFs.comb$Enh.class=="NPC_prim"]])
writeLines(human.lookup$Gene.name[human.lookup$Gene.stable.ID%in%h.TFs.comb$Gene.ID[h.TFs.comb$Enh.class=="7W_brain_prim"]])

writeLines(m.lookup$Gene.name[m.lookup$Gene.stable.ID%in%ect.TFs$Gene.ID])


if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("VennDetail")
library(VennDetail)

mECT=toupper(m.lookup$Gene.name[m.lookup$Gene.stable.ID%in%ect.TFs$Gene.ID])
NPC=toupper(human.lookup$Gene.name[human.lookup$Gene.stable.ID%in%h.TFs.comb$Gene.ID[h.TFs.comb$Enh.class=="NPC_prim"]])
h.brain=toupper((human.lookup$Gene.name[human.lookup$Gene.stable.ID%in%h.TFs.comb$Gene.ID[h.TFs.comb$Enh.class=="7W_brain_prim"]]))

plot(venndetail(list(mECT=mECT,NPC=NPC,h.Brain=h.brain)))
mECT
NPC
dev.off()


####Human####

setwd("G:/Headstone/Enh_priming/")

library(ggplot2)
library(tidyverse)
library(qlcMatrix)
library(ggpubr)

TFs=read.delim("./Annotations/Human/TFs.txt",h=F)[,1]
Gene.names=read.delim("./Annotations/Human/ENSEMBL_IDs_to_gene_names.txt")
head(Gene.names)



##Loading protein coding genes

tss=read.delim("./Annotations/Human/ENSEMBL_coding_gene_TSS.txt",h=F)
gene.ids=unique(tss[,4])

#Getting TF ENSEMBL IDs
TFs
Gene.names=read.delim("./Annotations/Human/ENSEMBL_IDs_to_gene_names.txt")
h.TF.genes=Gene.names[Gene.names$Gene.name%in%TFs,]
head(h.TF.genes)

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



#Getting upregulated genes in NPC transition
head(xie)
hist(xie$NPC)
xie$npc.log2fc=xie$NPC-xie$H1
expressed.npc=xie$ID[xie$NPC>0]
expressed.msc=xie$ID[xie$MSC>0]
expressed.me=xie$ID[xie$ME>0]


head(expressed.npc)
npc.TFs=Gene.names$Gene.name[Gene.names$Gene.stable.ID%in%comb$Gene.ID[comb$Enh.class=="NPC_prim"&comb$Gene.ID%in%expressed.npc]]
msc.TFs=Gene.names$Gene.name[Gene.names$Gene.stable.ID%in%comb$Gene.ID[comb$Enh.class=="MSC_prim"&comb$Gene.ID%in%expressed.msc]]
me.TFs=Gene.names$Gene.name[Gene.names$Gene.stable.ID%in%comb$Gene.ID[comb$Enh.class=="ME_prim"&comb$Gene.ID%in%expressed.me]]
write.table(npc.TFs,"./TF_network_analysis/hNPC_prim_TFs.txt",sep="\n",col.names = F,row.names = F,quote = F)
write.table(msc.TFs,"./TF_network_analysis/hMSC_prim_TFs.txt",sep="\n",col.names = F,row.names = F,quote = F)
write.table(me.TFs,"./TF_network_analysis/hME_prim_TFs.txt",sep="\n",col.names = F,row.names = F,quote = F)


xie$npc.log2fc=xie$NPC-xie$H1
expressed.npc=xie$ID[xie$NPC>0]



head(comb)

npc.upregulated=xie$ID[xie$npc.log2fc>2]
npc.TFs=h.TF.genes$Gene.stable.ID[h.TF.genes$Gene.stable.ID%in%expressed.npc]
up.npc.TFs=h.TF.genes$Gene.stable.ID[h.TF.genes$Gene.stable.ID%in%expressed.npc&h.TF.genes$Gene.stable.ID%in%npc.upregulated]
npc.enh.genes=comb$Gene.ID[comb$Enh.class=="NPC_prim"]
up.npc.enh.gens=h.TF.genes$Gene.name[h.TF.genes$Gene.stable.ID%in%(npc.enh.genes[npc.enh.genes%in%up.npc.TFs])]
writeLines(up.npc.enh.gens)

write.table(up.npc.enh.gens,"./TF_network_analysis/hNPC_prim_upTFs.txt",sep="\n",col.names = F,row.names = F,quote = F)

xie$ME.log2fc=xie$ME-xie$H1
expressed.me=xie$ID[xie$ME>0]
me.upregulated=xie$ID[xie$ME.log2fc>2]
up.me.TFs=h.TF.genes$Gene.stable.ID[h.TF.genes$Gene.stable.ID%in%expressed.me&h.TF.genes$Gene.stable.ID%in%me.upregulated]
me.enh.genes=comb$Gene.ID[comb$Enh.class=="ME_prim"]
up.me.enh.gens=h.TF.genes$Gene.name[h.TF.genes$Gene.stable.ID%in%(me.enh.genes[me.enh.genes%in%up.me.TFs])]
writeLines(up.me.enh.gens)
write.table(up.me.enh.gens,"./TF_network_analysis/hME_prim_upTFs.txt",sep="\n",col.names = F,row.names = F,quote = F)



xie$MSC.log2fc=xie$MSC-xie$H1
expressed.MSC=xie$ID[xie$MSC>0]
MSC.upregulated=xie$ID[xie$MSC.log2fc>2]
up.MSC.TFs=h.TF.genes$Gene.stable.ID[h.TF.genes$Gene.stable.ID%in%expressed.MSC&h.TF.genes$Gene.stable.ID%in%MSC.upregulated]
MSC.enh.genes=comb$Gene.ID[comb$Enh.class=="MSC_prim"]
up.MSC.enh.gens=h.TF.genes$Gene.name[h.TF.genes$Gene.stable.ID%in%(MSC.enh.genes[MSC.enh.genes%in%up.MSC.TFs])]
writeLines(up.MSC.enh.gens)
write.table(up.MSC.enh.gens,"./TF_network_analysis/hMSC_prim_upTFs.txt",sep="\n",col.names = F,row.names = F,quote = F)


expressed.cerb=data.coding$ID[data.coding$Tissue=="Fetal_Cerebrum"&data.coding$percentile>44]
fbrain.tfs=h.TF.genes$Gene.stable.ID[h.TF.genes$Gene.stable.ID%in%expressed.cerb]
fbrain.genes=comb$Gene.ID[comb$Enh.class=="7W_brain_prim"]
expr.fbrain.genes=h.TF.genes$Gene.name[h.TF.genes$Gene.stable.ID%in%(fbrain.genes[fbrain.genes%in%fbrain.tfs])]
writeLines(expr.fbrain.genes)
write.table(expr.fbrain.genes,"./TF_network_analysis/hBrain_prim_TFs.txt",sep="\n",col.names = F,row.names = F,quote = F)


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



####Mouse####


setwd("E:/Headstone/Enh_priming/")
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
  
  
  
  fetal.comb=merge(fetal,fetal.heart,by="Gene.ID")
  head(fetal.comb)  
  hist(fetal.comb$E11_Brain.fore)
  hist(fetal.comb$E11_Heart)
  fetal.comb$fc.brain=fetal.comb$E11_Brain.fore-fetal.comb$E11_Liver
  
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
}

expr.scnmt=scnmt[expr.scnmt,]
head(expr.scnmt)
hist(expr.scnmt$E7.5_Ectoderm)
expr.scnmt$ect.fc=expr.scnmt$E7.5_Ectoderm-expr.scnmt$E4.5_Epiblast
expr.scnmt$end.fc=expr.scnmt$E7.5_Endoderm-expr.scnmt$E4.5_Epiblast
expr.scnmt$mes.fc=expr.scnmt$E7.5_Mesoderm-expr.scnmt$E4.5_Epiblast

hist(expr.scnmt$ect.fc)
up.ect.genes=expr.scnmt$Gene.ID[expr.scnmt$ect.fc>0]
up.mes.genes=expr.scnmt$Gene.ID[expr.scnmt$mes.fc>0]
up.end.genes=expr.scnmt$Gene.ID[expr.scnmt$end.fc>0]

up.brain.genes=fetal.comb$Gene.ID[fetal.comb$E11_Brain.fore>1&fetal.comb$fc.brain>0]


TFs=read.delim("./Annotations/Human/TFs.txt",h=F)
id.lookup=read.delim("./Annotations/Mouse/ENSEMBL_IDs_to_gene_names.txt")
id.lookup$Gene.name=toupper(id.lookup$Gene.name)
tf.lookup=id.lookup[id.lookup$Gene.name%in%TFs$V1,]

ect.enh.tfs=tf.lookup$Gene.name[tf.lookup$Gene.stable.ID%in%comb$Gene.ID[comb$Enh.class=="mECT_prim"]]
mes.enh.tfs=tf.lookup$Gene.name[tf.lookup$Gene.stable.ID%in%comb$Gene.ID[comb$Enh.class=="mMES_prim"]]
end.enh.tfs=tf.lookup$Gene.name[tf.lookup$Gene.stable.ID%in%comb$Gene.ID[comb$Enh.class=="mEND_prim"]]
brain.enh.tfs=tf.lookup$Gene.name[tf.lookup$Gene.stable.ID%in%comb$Gene.ID[comb$Enh.class=="brain_prim"]]

up.ect.tf=ect.enh.tfs[ect.enh.tfs%in%tf.lookup$Gene.name[tf.lookup$Gene.stable.ID%in%up.ect.genes]]
up.mes.tf=mes.enh.tfs[mes.enh.tfs%in%tf.lookup$Gene.name[tf.lookup$Gene.stable.ID%in%up.mes.genes]]
up.end.tf=end.enh.tfs[end.enh.tfs%in%tf.lookup$Gene.name[tf.lookup$Gene.stable.ID%in%up.end.genes]]
up.brain.tf=brain.enh.tfs[brain.enh.tfs%in%tf.lookup$Gene.name[tf.lookup$Gene.stable.ID%in%up.brain.genes]]

write.table(ect.enh.tfs,"./TF_network_analysis/mECT_prim_TFs.txt",sep="\n",col.names = F,row.names = F,quote = F)
write.table(mes.enh.tfs,"./TF_network_analysis/mMES_prim_TFs.txt",sep="\n",col.names = F,row.names = F,quote = F)
write.table(end.enh.tfs,"./TF_network_analysis/mEND_prim_TFs.txt",sep="\n",col.names = F,row.names = F,quote = F)
write.table(brain.enh.tfs,"./TF_network_analysis/mBrain_prim_TFs.txt",sep="\n",col.names = F,row.names = F,quote = F)


write.table(up.ect.tf,"./TF_network_analysis/mECT_prim_upTFs.txt",sep="\n",col.names = F,row.names = F,quote = F)
write.table(up.mes.tf,"./TF_network_analysis/mMES_prim_upTFs.txt",sep="\n",col.names = F,row.names = F,quote = F)
write.table(up.end.tf,"./TF_network_analysis/mEND_prim_upTFs.txt",sep="\n",col.names = F,row.names = F,quote = F)
write.table(up.brain.tf,"./TF_network_analysis/mBrain_prim_upTFs.txt",sep="\n",col.names = F,row.names = F,quote = F)




####Plotting TF numbers####
mect.tf=read.delim("./TF_network_analysis/mECT_prim_TFs.txt",h=F)[,1]
mmes.tf=read.delim("./TF_network_analysis/mMES_prim_TFs.txt",h=F)[,1]
mend.tf=read.delim("./TF_network_analysis/mEND_prim_TFs.txt",h=F)[,1]
mBrain.tf=read.delim("./TF_network_analysis/mBrain_prim_TFs.txt",h=F)[,1]

mect.tf.up=read.delim("./TF_network_analysis/mECT_prim_upTFs.txt",h=F)[,1]
mmes.tf.up=read.delim("./TF_network_analysis/mMES_prim_upTFs.txt",h=F)[,1]
mend.tf.up=read.delim("./TF_network_analysis/mEND_prim_upTFs.txt",h=F)[,1]
mBrain.tf.up=read.delim("./TF_network_analysis/mBrain_prim_upTFs.txt",h=F)[,1]

hnpc.tf=read.delim("./TF_network_analysis/hNPC_prim_TFs.txt",h=F)[,1]
hme.tf=read.delim("./TF_network_analysis/hME_prim_TFs.txt",h=F)[,1]
hmsc.tf=read.delim("./TF_network_analysis/hMSC_prim_TFs.txt",h=F)[,1]
hbrain.tf=read.delim("./TF_network_analysis/hBrain_prim_TFs.txt",h=F)[,1]

hnpc.tf.up=read.delim("./TF_network_analysis/hNPC_prim_upTFs.txt",h=F)[,1]
hme.tf.up=read.delim("./TF_network_analysis/hME_prim_upTFs.txt",h=F)[,1]
hmsc.tf.up=read.delim("./TF_network_analysis/hMSC_prim_upTFs.txt",h=F)[,1]


library(VennDetail)

dev.off()
plot(venndetail(list(mECT=mect.tf,mMES=mmes.tf,mEND=mend.tf,mBrain=mBrain.tf)))

dev.off()
plot(venndetail(list(mECT=mect.tf.up,mMES=mmes.tf.up,mEND=mend.tf.up,mBrain=mBrain.tf.up)))

dev.off()
plot(venndetail(list(hNPC=hnpc.tf,hME=hme.tf,hMSC=hmsc.tf,hBrain=hbrain.tf)))
dev.off()
plot(venndetail(list(hNPC=hnpc.tf.up,hME=hme.tf.up,hMSC=hmsc.tf.up,hBrain=hbrain.tf)))

dev.off()
plot(venndetail(list(hNPC=hnpc.tf.up,mECT=mect.tf,mBrain=mBrain.tf,hBrain=hbrain.tf)))

dev.off()
plot(venndetail(list(mMES=mmes.tf,mEND=mend.tf,hME=hme.tf,hMSC=hmsc.tf)))


writeLines(mect.tf.up)
