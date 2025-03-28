setwd("G:/Headstone/Enh_priming/")


m.enh=read.delim("./Annotations/Mouse/Comb_Enh_class_lookup.txt",h=F)
head(m.enh)

m.enh.coord=read.delim("./Annotations/Mouse/Enhancers/All_tissue_enh.txt",h=F)
head(m.enh.coord)


summary(as.factor(m.enh.coord$V5))
summary(as.factor(m.enh$V2))

summary(m.enh$V1%in%m.enh.coord$V4)
m.enh[!m.enh$V1%in%m.enh.coord$V4,]

class.of.interest=c("forebrain_prim","forebrain_nonprim",
                    "mECT_prim","mECT_nonprim",
                    "mEND_prim","mEND_nonprim",
                    "mMES_prim","mMES_nonprim",
                    "heart_prim","heart_nonprim",
                    "liver_prim","liver_nonprim")

m.enh.sub=m.enh[m.enh$V2%in%class.of.interest,]
colnames(m.enh.sub)=c("Enh.ID","Enh.class")

colnames(m.enh.coord)=c("chr","start","end","Enh.ID","tissue")
head(m.enh.coord)
tog=merge(m.enh.coord,m.enh.sub,by="Enh.ID")
head(tog)
tog=tog[,c(2:4,1,6)]



####Human####
setwd("G:/Headstone/Enh_priming/")
library(ggplot2)
library(ggalluvial)

#Getting All human enhancer groups (before lineage specific subsetting)

xie.enh=read.delim("./Annotations/Human/Xie_enh_coord_hg38.txt",h=F)
xie.enh[,5]=gsub("_all_enhancers","",xie.enh[,5])
xie.enh[,5]=gsub("_enriched_enhancers","",xie.enh[,5])
brain.enh=read.delim("./Annotations/Human/Enhancers/7W_Brain_Enh.txt",h=F)
brain.enh[,5]=rep("brain")
fetal.enh=read.delim("./Annotations/Human/Enhancers/Fetal_enhancers.txt",h=F)
human.enh=do.call(rbind,list(xie.enh,brain.enh,fetal.enh))
human.enh=unique(human.enh[,4:5])
colnames(human.enh)=c("ID","Enh.class")

classes=unique(human.enh$Enh.class)

comb.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(comb.lookup)=c("ID","Enh.class")

#loading the RPKM values
xie.rpkms=read.delim("./Data/human/ChIP/Xie_Enh_1_5kb_ChIP_logRPKM_v2.txt")[,c(1,13:18)]
fetal.rpkms=read.delim("./Data/human/ChIP/Fetal_Enh_1_5kb_ChIP_logRPKM.txt")[,c(1,13:18)]
brain.rpkms=read.delim("./Data/human/ChIP/7W_Brain_Enh_1_5kb_ChIP_logRPKM_v2.txt")[,c(1,13:18)]
tog.rpkms=do.call(rbind,list(xie.rpkms,fetal.rpkms,brain.rpkms))

set.k4me1=1.2
set.k27ac=2.97
set.k27me3=2

annotation.files=list.files("./Annotations/Human/Enhancers/",pattern = "enh.txt")


poised=tog.rpkms$Probe[tog.rpkms$H1_K27ac<set.k27ac&tog.rpkms$H1_K27me3>=set.k27me3]
poised.ids=human.enh[human.enh$ID%in%poised,]
all.enh=do.call(rbind,list(xie.enh,brain.enh,fetal.enh))
poised.enh=all.enh[all.enh$V4%in%poised.ids$ID,]
poised.enh$V5=paste0(poised.enh$V5,"_poised")

write.table(poised.enh,"./Rebuttal/Human_Poised_enhancers.txt",sep="\t",col.names = F,row.names = F,quote = F)
specific.poised=poised.enh[poised.enh$V4%in%comb.lookup$ID,]
write.table(specific.poised,"./Rebuttal/Human_Specific_Poised_enhancers.txt",sep="\t",col.names = F,row.names = F,quote = F)


####Mouse####

setwd("G:/Headstone/Enh_priming/")
library(ggplot2)
library(ggalluvial)



set.k4me1=-0.4
set.k27ac=1.2
set.k27me3=1.25

all.enh=read.delim("./Annotations/Mouse/Enhancers/All_tissue_enh.txt",h=F)

mouse.enh=all.enh[,4:5]
colnames(mouse.enh)=c("ID","Enh.class")
mouse.enh$Enh.class=gsub("Nonspecific_","",mouse.enh$Enh.class)

classes=unique(mouse.enh$Enh.class)

comb.lookup=read.delim("./Annotations/Mouse/Comb_Enh_class_lookup.txt",h=F)
colnames(comb.lookup)=c("ID","Enh.class")

tog.rpkms=read.delim("./Data/mouse/ChIP/All_Enh_Aled_TETKO_ChIP_logRPKM.txt")[,c(7,13:15)]
colnames(tog.rpkms)=c("ID","K27ac","K27me3","K4me1")

poised=tog.rpkms$ID[tog.rpkms$K27ac<3&tog.rpkms$K27me3>=set.k27me3]

poised.enh=all.enh[all.enh$V4%in%poised,]
poised.enh$V5=gsub("Nonspecific_","",poised.enh$V5)
poised.enh$V5=paste0(poised.enh$V5,"_poised")
write.table(poised.enh,"./Rebuttal/Mouse_Poised_enhancers.txt",sep="\t",col.names = F,row.names = F,quote = F)
specific.poised=poised.enh[poised.enh$V4%in%comb.lookup$ID,]
write.table(specific.poised,"./Rebuttal/Mouse_Specific_Poised_enhancers.txt",sep="\t",col.names = F,row.names = F,quote = F)
head(specific.poised)
hpoised



setwd("G:/Headstone/Enh_priming/")
library(ggplot2)
library(tidyverse)

hpoised=read.delim("./Rebuttal/Human_Poised_enhancers.txt",h=F)[,4:5]
colnames(hpoised)=c("ID","Enh.class")


hpoised=hpoised[hpoised$ID%in%z$V4,]
####Stat tests for Theunissen human primed DNA met (Figure 1)####

{to.plot=c("Global","hESC_specific","NPC_prim","NPC_nonprim","NPC_poised","ME_prim","ME_nonprim","ME_poised","7W_brain_prim","7W_brain_nonprim","brain_poised")
enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(enh.lookup)=c("ID","Enh.class")
enh.lookup=rbind(enh.lookup,hpoised)
hesc=read.delim("./Data/human/DNA_methyl/Theunissen_2016_hESC_WGBS.txt")[,c(7,13:14)]
poised.hesc=read.delim("./Rebuttal/Human_poised_theunissen_hesc_WGBS.txt")[,c(1,13:14)]
colnames(poised.hesc)=colnames(hesc)
hesc=rbind(hesc,poised.hesc)

hesc=merge(hesc,enh.lookup,by="ID")

global.hesc=read.delim("./Data/human/DNA_methyl/Theunissen_2016_hESC_WGBS_Global.txt")[,c(1,13:14)]
global.hesc$Enh.class=rep("Global")
colnames(global.hesc)=colnames(hesc)
hesc=rbind(hesc,global.hesc)
hesc=pivot_longer(hesc[,c(4,2:3)],cols = -1,names_to = "Tissue",values_to = "Methyl")
hesc=hesc[hesc$Enh.class%in%to.plot,]
hesc$Enh.class=factor(hesc$Enh.class,levels = to.plot)

sub.df=hesc[hesc$Tissue=="Primed_hESC",]
}
summary(aov(Methyl ~ Enh.class,sub.df[sub.df$Enh.class%in%c("NPC_prim","NPC_nonprim"),]))
summary(aov(Methyl ~ Enh.class,sub.df[sub.df$Enh.class%in%c("ME_prim","ME_nonprim"),]))
rm(hesc,sub.df,global.hesc,enh.lookup,to.plot)
library(ggplot2)
head(sub.df)
ggplot(sub.df,aes(x=Enh.class,y=Methyl,fill=Enh.class))+geom_boxplot()


####Human Homer Heatmaps####

setwd("G:/Headstone/Enh_priming/")

library(ComplexHeatmap)
library(circlize)
library("progress")
library(ggplot2)

lookup=read.delim("./Rebuttal/Human_Specific_Poised_enhancers.txt",h=F)[,4:5]
#lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(lookup)=c("ID","class")
lookup=lookup[!is.na(lookup$ID),]
classes=unique(lookup$class)

hist.dir="./Homer/Homer_hist/human/"
hist.files=list.files(hist.dir)

hist.files
colour.list=c("H1_K27ac"="green4",
              "H1_K27me3"="green4",
              "H1_K4me1"="green4",
              "ME_K27ac"="blue",
              "ME_K4me1"="blue",
              "MSC_K27ac"="red",
              "MSC_K4me1"="red",
              "NPC_K27ac"="purple4",
              "NPC_K4me1"="purple4")

hist.list=c("H1_K27ac"="H1_K27ac_hist.txt",
            "H1_K27me3"="H1_K27me3_hist.txt",
            "H1_K4me1"="H1_K4me1_hist.txt",
            "ME_K27ac"="ME_K27ac_hist.txt",
            "ME_K4me1"="ME_K4me1_hist.txt",
            "MSC_K27ac"="MSC_K27ac_hist.txt",
            "MSC_K4me1"="MSC_K4me1_hist.txt",
            "NPC_K27ac"="NPC_K27ac_hist.txt",
            "NPC_K4me1"="NPC_K4me1_hist.txt")




ref.hist.file=hist.files[grepl("H1_K4me1",hist.files)]
ref.hist=read.delim(paste0(hist.dir,ref.hist.file))
ref.mat=as.matrix(ref.hist[,2:ncol(ref.hist)])
row.names(ref.mat)=ref.hist[,1]
not.in.ref=unique(as.character(lookup$ID[lookup$ID%in%rownames(ref.mat)==F]))
not.in.ref.mat=matrix(0,nrow=length(not.in.ref),ncol = ncol(ref.mat))
row.names(not.in.ref.mat)=not.in.ref
ref.mat=rbind(ref.mat,not.in.ref.mat)
ref.mat=ref.mat[order(rowSums(ref.mat[,40:60]),decreasing = T),]
ref.ids=rownames(ref.mat)[rownames(ref.mat)%in%lookup$ID]

rm(ref.mat,ref.hist,not.in.ref.mat,ref.hist.file,not.in.ref)



saturation.list=c("H1_K27ac"=0.08,
                  "H1_K27me3"=0.1,
                  "H1_K4me1"=0.08,
                  "ME_K27ac"=0.08,
                  "ME_K4me1"=0.08,
                  "MSC_K27ac"=0.08,
                  "MSC_K4me1"=0.08,
                  "NPC_K27ac"=0.08,
                  "NPC_K4me1"=0.08)

pb <- progress_bar$new(
  format = "  Getting heatmaps [:bar] :percent eta: :eta",
  total = (length(classes)*length(hist.files)), clear = FALSE, width= 60)




for(i in names(hist.list)){
  hist=read.delim(paste0(hist.dir,hist.list[i]))
  mat=as.matrix(hist[,2:ncol(hist)])
  row.names(mat)=hist[,1]
  not.in.hist=unique(as.character(lookup$ID[lookup$ID%in%rownames(mat)==F]))
  not.mat=matrix(0,ncol = ncol(mat),nrow=length(not.in.hist))
  row.names(not.mat)=not.in.hist
  mat=rbind(mat,not.mat)
  mat=mat[ref.ids,]
  col.ramp=colorRamp2(c(min(mat),saturation.list[i]),c("white",colour.list[i]))
  col.ramp2=colorRamp2(c(min(mat),saturation.list[i]),c("blue","red"))
  for(n in classes){
    ids=lookup$ID[lookup$class==n]
    mat2=mat[rownames(mat)%in%ids,]
    hm.plot=Heatmap(mat2,cluster_rows=F,cluster_columns=F,
                    show_row_names=F,show_column_names=F,
                    show_heatmap_legend=F,col=col.ramp,use_raster = T)
    hm.plot2=Heatmap(mat2,cluster_rows=F,cluster_columns=F,
                     show_row_names=F,show_column_names=F,
                     show_heatmap_legend=F,col=col.ramp2,use_raster = T)
    output.name=paste0("./Plots/Heatmaps/human/",i,"_",n,".jpg")
    jpeg(output.name,w=1.5,h=5,units='in',res=600)
    draw(hm.plot)
    dev.off()
    output.name=paste0("./Plots/Heatmaps/human/RB/",i,"_",n,".jpg")
    jpeg(output.name,w=1.5,h=5,units='in',res=600)
    draw(hm.plot2)
    dev.off()
    pb$tick()
  }
}

####Plotting mouse heatmaps####
setwd("G:/Headstone/Enh_priming/")

library(ComplexHeatmap)
library(circlize)
library("progress")
library(ggplot2)

lookup=read.delim("./Rebuttal/Mouse_Specific_Poised_enhancers.txt",h=F)[,4:5]
colnames(lookup)=c("ID","class")
lookup=lookup[!is.na(lookup$ID),]
classes=unique(lookup$class)

hist.dir="./Homer/Homer_hist/mouse/"
hist.files=list.files(hist.dir)


colour.list=c(
  "Aled_2i_H3K27ac"="green4",  "Aled_2i_H3K27me3"="green4",  "Aled_2i_H3K4me1"="green4",
  "Aled_EpiLC_H3K27ac"="green4",  "Aled_EpiLC_H3K27me3"="green4",  "Aled_EpiLC_H3K4me1"="green4",
  "E11_forebrian_H3K27ac"="black",
  "E11_heart_H3K27ac"="black",
  "E11_hindbrain_H3K27ac"="black",
  "E11_liver_H3K27ac"="black",
  "E11_midbrain_H3K27ac"="black",
  "mECT_H3K27ac"="purple4",  "mEND_H3K27ac"="blue",  "mMES_H3K27ac"="red"
)
hist.list=c(
  "Aled_2i_H3K27ac"="Aled_2i_H3K27ac_hist.txt",  "Aled_2i_H3K27me3"="Aled_2i_H3K27me3_hist.txt",  "Aled_2i_H3K4me1"="Aled_2i_H3K4me1_hist.txt",
  "Aled_EpiLC_H3K27ac"="Aled_EpiLC_H3K27ac_hist.txt",  "Aled_EpiLC_H3K27me3"="Aled_EpiLC_H3K27me3_hist.txt",   "Aled_EpiLC_H3K4me1"="Aled_EpiLC_H3K4me1_hist.txt",
  "E11_forebrian_H3K27ac"="E11_forebrain_K27ac_hist.txt",
  "E11_heart_H3K27ac"="E11_heart_K27ac_hist.txt",
  "E11_hindbrain_H3K27ac"="E11_hindbrain_K27ac_hist.txt",
  "E11_liver_H3K27ac"="E11_liver_K27ac_hist.txt",
  "E11_midbrain_H3K27ac"="E11_midbrain_K27ac_hist.txt",
  "mECT_H3K27ac"="mECT_K27ac_hist.txt",  "mEND_H3K27ac"="mEND_K27ac_hist.txt",  "mMES_H3K27ac"="mMES_K27ac_hist.txt"
)



ref.hist.file=hist.files[grepl("2i_H3K4me1",hist.files)]
ref.hist=read.delim(paste0(hist.dir,ref.hist.file))
ref.mat=as.matrix(ref.hist[,2:ncol(ref.hist)])
row.names(ref.mat)=ref.hist[,1]
not.in.ref=unique(as.character(lookup$ID[lookup$ID%in%rownames(ref.mat)==F]))
not.in.ref.mat=matrix(0,nrow=length(not.in.ref),ncol = ncol(ref.mat))
row.names(not.in.ref.mat)=not.in.ref
ref.mat=rbind(ref.mat,not.in.ref.mat)
ref.mat=ref.mat[order(rowSums(ref.mat[,40:60]),decreasing = T),]
ref.ids=rownames(ref.mat)[rownames(ref.mat)%in%lookup$ID]
rm(ref.mat,ref.hist,not.in.ref.mat,ref.hist.file,not.in.ref)


get.heatmap=function(ids,mark,to.ref=T){
  hist=read.delim(paste0(hist.dir,hist.list[mark]))
  mat=as.matrix(hist[,2:ncol(hist)])
  row.names(mat)=hist[,1]
  not.in.hist=unique(as.character(lookup$ID[lookup$ID%in%rownames(mat)==F]))
  not.mat=matrix(0,ncol = ncol(mat),nrow=length(not.in.hist))
  row.names(not.mat)=not.in.hist
  mat=rbind(mat,not.mat)
  if(to.ref){
    mat=mat[ref.ids,]
  }else{
    mat=mat[order(rowSums(mat[,40:60]),decreasing=T),]
  }
  
  col.ramp=colorRamp2(c(min(mat),saturation.list[mark]),c("white",colour.list[mark]))
  mat=mat[rownames(mat)%in%ids,]
  hm.plot=Heatmap(mat,cluster_rows=F,cluster_columns=F,
                  show_row_names=F,show_column_names=F,
                  show_heatmap_legend=F,col=col.ramp,use_raster = T)
  return(hm.plot)
}
sat.test.heatmap=function(ids,mark,sat,to.ref=T){
  hist=read.delim(paste0(hist.dir,hist.list[mark]))
  mat=as.matrix(hist[,2:ncol(hist)])
  row.names(mat)=hist[,1]
  not.in.hist=unique(as.character(lookup$ID[lookup$ID%in%rownames(mat)==F]))
  not.mat=matrix(0,ncol = ncol(mat),nrow=length(not.in.hist))
  row.names(not.mat)=not.in.hist
  mat=rbind(mat,not.mat)
  if(to.ref){
    mat=mat[ref.ids,]
  }else{
    mat=mat[order(rowSums(mat[,40:60]),decreasing=T),]
  }
  
  col.ramp=colorRamp2(c(min(mat),sat),c("white",colour.list[mark]))
  mat=mat[rownames(mat)%in%ids,]
  hm.plot=Heatmap(mat,cluster_rows=F,cluster_columns=F,
                  show_row_names=F,show_column_names=F,
                  show_heatmap_legend=F,col=col.ramp,use_raster = T)
  return(hm.plot)
}
saturation.list=c(
  "Aled_2i_H3K27ac"=0.08,  "Aled_2i_H3K27me3"=0.08,  "Aled_2i_H3K4me1"=0.07,
  "Aled_EpiLC_H3K27ac"=0.08,  "Aled_EpiLC_H3K27me3"=0.08,  "Aled_EpiLC_H3K4me1"=0.07,
  "E11_forebrian_H3K27ac"=0.1,
  "E11_heart_H3K27ac"=0.1,
  "E11_hindbrain_H3K27ac"=0.1,
  "E11_liver_H3K27ac"=0.1,
  "E11_midbrain_H3K27ac"=0.1,
  "mECT_H3K27ac"=0.1,  "mEND_H3K27ac"=0.1,  "mMES_H3K27ac"=0.1)

pb <- progress_bar$new(
  format = "  Getting heatmaps [:bar] :percent eta: :eta",
  total = (length(classes)*length(hist.files)), clear = FALSE, width= 60)

for(i in names(hist.list)){
  hist=read.delim(paste0(hist.dir,hist.list[i]))
  mat=as.matrix(hist[,2:ncol(hist)])
  row.names(mat)=hist[,1]
  not.in.hist=unique(as.character(lookup$ID[lookup$ID%in%rownames(mat)==F]))
  not.mat=matrix(0,ncol = ncol(mat),nrow=length(not.in.hist))
  row.names(not.mat)=not.in.hist
  mat=rbind(mat,not.mat)
  mat=mat[ref.ids,]
  col.ramp=colorRamp2(c(min(mat),saturation.list[i]),c("white",colour.list[i]))
  col.ramp2=colorRamp2(c(min(mat),saturation.list[i]),c("blue","red"))
  for(n in classes){
    ids=lookup$ID[lookup$class==n]
    mat2=mat[rownames(mat)%in%ids,]
    hm.plot=Heatmap(mat2,cluster_rows=F,cluster_columns=F,
                    show_row_names=F,show_column_names=F,
                    show_heatmap_legend=F,col=col.ramp,use_raster = T)
    hm.plot2=Heatmap(mat2,cluster_rows=F,cluster_columns=F,
                     show_row_names=F,show_column_names=F,
                     show_heatmap_legend=F,col=col.ramp2,use_raster = T)
    output.name=paste0("./Plots/Heatmaps/mouse/",i,"_",n,".jpg")
    jpeg(output.name,w=1.5,h=5,units='in',res=600)
    draw(hm.plot)
    dev.off()
    output.name=paste0("./Plots/Heatmaps/mouse/RB/",i,"_",n,".jpg")
    jpeg(output.name,w=1.5,h=5,units='in',res=600)
    draw(hm.plot2)
    dev.off()
    pb$tick()
  }
}
