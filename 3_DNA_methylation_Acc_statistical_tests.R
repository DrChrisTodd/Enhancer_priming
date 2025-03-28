setwd("G:/Headstone/Enh_priming/")
library(ggplot2)
library(tidyverse)



####Stat tests for Theunissen human primed DNA met (Figure 1)####

{to.plot=c("Global","hESC_specific","NPC_prim","NPC_nonprim","ME_prim","ME_nonprim")
enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(enh.lookup)=c("ID","Enh.class")
hesc=read.delim("./Data/human/DNA_methyl/Theunissen_2016_hESC_WGBS.txt")[,c(7,13:14)]
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

####stat test for Shirane mouse ESC DNA met (Figure 1)####
{mouse.lookup=read.delim("./Annotations/Mouse/Comb_Enh_class_lookup.txt",h=F)
colnames(mouse.lookup)=c("ID","Enh.class")
to.plot=c("Global","mESC_enh","mECT_prim","mECT_nonprim","mMES_prim","mMES_nonprim","mEND_prim","mEND_nonprim")
mesc=read.delim("./Data/mouse/DNA_methyl/Shirane_2016_WGBS.txt")[,c(7,14:15)]
mesc.global=read.delim("./Data/mouse/DNA_methyl/Shirane_2016_WGBS_Global.txt")[,c(1,14:15)]
colnames(mesc.global)=colnames(mesc)
mesc=merge(mesc,mouse.lookup,by="ID")
mesc.global$Enh.class=rep("Global")
mesc=rbind(mesc,mesc.global)
mesc=pivot_longer(mesc[,c(4,2:3)],cols=-1,names_to = "Tissue",values_to = "Methyl")
mesc=mesc[mesc$Enh.class%in%to.plot,]
mesc$Enh.class=factor(mesc$Enh.class,levels=to.plot)
}
summary(aov(Methyl ~ Enh.class,mesc[mesc$Enh.class%in%c("mECT_prim","mECT_nonprim")&mesc$Tissue=="EpiLC",]))
summary(aov(Methyl ~ Enh.class,mesc[mesc$Enh.class%in%c("mMES_prim","mMES_nonprim")&mesc$Tissue=="EpiLC",]))

summary(aov(Methyl ~ Enh.class,mesc[mesc$Enh.class%in%c("mECT_prim","mECT_nonprim")&mesc$Tissue=="ESC",]))
summary(aov(Methyl ~ Enh.class,mesc[mesc$Enh.class%in%c("mMES_prim","mMES_nonprim")&mesc$Tissue=="ESC",]))


summary(aov(Methyl ~ Enh.class,mesc[mesc$Enh.class%in%c("mEND_prim","mEND_nonprim")&mesc$Tissue=="EpiLC",]))
summary(aov(Methyl ~ Enh.class,mesc[mesc$Enh.class%in%c("mEND_prim","mEND_nonprim")&mesc$Tissue=="ESC",]))


 rm(mesc,mesc.global,mouse.lookup)


####Stat tests for human embryonic scNMT seq data####
{
enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(enh.lookup)=c("ID","Enh.class")
hemb=read.delim("G:/Headstone/single_cell/human scNMT/metacc/boxplots_enhancers/boxplot_met_nich.txt")
}
summary(aov(rate~anno,hemb[hemb$lineage%in%c("ICM")&hemb$anno%in%c("NPC_prim","NPC_nonprim"),]))
summary(aov(rate~anno,hemb[hemb$lineage%in%c("Epiblast")&hemb$anno%in%c("NPC_prim","NPC_nonprim"),]))
summary(aov(rate~anno,hemb[hemb$lineage%in%c("ICM")&hemb$anno%in%c("ME_prim","ME_nonprim"),]))
summary(aov(rate~anno,hemb[hemb$lineage%in%c("Epiblast")&hemb$anno%in%c("ME_prim","ME_nonprim"),]))
rm(enh.lookup,hemb)

####Stat test for mouse embryonic scNMT data####
memb=read.delim("G:/Headstone/single_cell/scnmt_gastrulation-develop/metacc/boxplots_enhancers/met_values.txt")
summary(aov(rate~anno,memb[memb$stage_lineage=="E4.5 Epiblast"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="E4.5 Epiblast"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="E6.5 Epiblast"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="E6.5 Epiblast"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="E4.5 Epiblast"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="E6.5 Epiblast"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

 rm(memb)

####Stat test for mouse embryonic COOLseq data####
memb=read.delim("G:/Headstone/single_cell/COOLseq/metacc/boxplots_enhancers/met_values_mouse.txt")[,c(1,2,3,4)]
global=read.delim("G:/Headstone/single_cell/COOLseq/metacc/boxplots_enhancers/met_global_values.txt")
colnames(global)[4]="rate"
global=global[,colnames(memb)]
memb=rbind(memb,global)
summary(aov(rate~anno,memb[memb$stage_lineage=="Zygote"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Zygote"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Zygote"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Zygote"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="2-cell"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="2-cell"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="2-cell"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="2-cell"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="4-cell"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="4-cell"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="4-cell"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="4-cell"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="8-cell"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="8-cell"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="8-cell"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="8-cell"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="Morula"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Morula"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Morula"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Morula"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="ICM"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="ICM"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="ICM"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="ICM"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

 summary(aov(rate~anno,memb[memb$stage_lineage=="TE"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="TE"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="TE"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="TE"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

memb=read.delim("G:/Headstone/single_cell/COOLseq/metacc/boxplots_enhancers/met_values_encode.txt")
summary(as.factor(memb$anno))

summary(aov(rate~anno,memb[memb$stage_lineage=="Zygote"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="2-cell"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="4-cell"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="8-cell"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Morula"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="ICM"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
 summary(aov(rate~anno,memb[memb$stage_lineage=="TE"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))

rm(memb)

####Stat test for mouse embryonic COOLseq Acc data####
memb=read.delim("G:/Headstone/single_cell/COOLseq/metacc/boxplots_enhancers/acc_values_mouse.txt")[,c(1:4)]
global=read.delim("G:/Headstone/single_cell/COOLseq/metacc/boxplots_enhancers/acc_global_values.txt")
colnames(global)[4]="rate"
global=global[,colnames(memb)]
memb2=read.delim("G:/Headstone/single_cell/COOLseq/metacc/boxplots_enhancers/acc_values_encode.txt")[,c(1:4)]
memb=rbind(memb,memb2)
memb=rbind(memb,global)
head(memb)

summary(aov(rate~anno,memb[memb$stage_lineage=="Zygote"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Zygote"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Zygote"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Zygote"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="2-cell"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="2-cell"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="2-cell"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="2-cell"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="4-cell"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="4-cell"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="4-cell"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="4-cell"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="8-cell"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="8-cell"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="8-cell"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="8-cell"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

 summary(aov(rate~anno,memb[memb$stage_lineage=="Morula"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Morula"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Morula"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Morula"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="ICM"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="ICM"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="ICM"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="ICM"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="TE"&memb$anno%in%c("all","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="TE"&memb$anno%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="TE"&memb$anno%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="TE"&memb$anno%in%c("mEND_prim","mEND_nonprim"),]))

summary(aov(rate~anno,memb[memb$stage_lineage=="Zygote"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
 summary(aov(rate~anno,memb[memb$stage_lineage=="2-cell"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="4-cell"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="8-cell"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="Morula"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="ICM"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))
summary(aov(rate~anno,memb[memb$stage_lineage=="TE"&memb$anno%in%c("forebrain_prim","forebrain_nonprim"),]))


####Stat test for human TET and DNMT KO in ESC methyl data####
setwd("G:/Headstone/Enh_priming/")
enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(enh.lookup)=c("ID","Enh.class")

to.plot=c("hESC_specific","MSC_specific","NPC_prim","NPC_nonprim","ME_prim","ME_nonprim")
df=read.delim("./Data/human/DNA_methyl/Charlston_2020_DNMT_TET_KO_WGBS.txt")
df=df[,c(7,13:16,20:22)]
df=merge(df,enh.lookup,by="ID")
plot.order=c("HUES8_WT","HUES8_DKO_P6","HUES8_TKO","HUES8_PKO","HUES8_PKO_TET1s_rescue","HUES8_PKO_TET2_rescue","HUES8_PKO_TET3_rescue")
df=pivot_longer(df[,c(9,2:8)],cols=-1,names_to="Tissue",values_to="Methyl")
df$Tissue=factor(df$Tissue,levels=plot.order)


TukeyHSD(aov(Methyl~Tissue,df[df$Enh.class=="NPC_nonprim",]))
 TukeyHSD(aov(Methyl~Tissue,df[df$Enh.class=="NPC_prim",]))
TukeyHSD(aov(Methyl~Tissue,df[df$Enh.class=="ME_nonprim",]))
 TukeyHSD(aov(Methyl~Tissue,df[df$Enh.class=="NPC_prim",]))
TukeyHSD(aov(Methyl~Tissue,df[df$Enh.class=="MSC_specific",]))

global=read.delim("./Data/human/DNA_methyl/Charlston_2020_DNMT_TET_KO_WGBS_Global.txt")
global=global[,c(7,13:16,20:22)]
global=pivot_longer(global,cols=-1,names_to = "Tissue",values_to="Methyl")
global$Tissue=factor(global$Tissue,levels=plot.order)
global$Enh.class=rep("Global")
global=global[,colnames(df)]
df2=rbind(df,global)

summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("HUES8_PKO_TET1s_rescue")&df2$Enh.class%in%c("ME_nonprim","Global"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("HUES8_PKO_TET1s_rescue")&df2$Enh.class%in%c("NPC_nonprim","Global"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("HUES8_PKO_TET2_rescue")&df2$Enh.class%in%c("NPC_nonprim","Global"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("HUES8_PKO_TET2_rescue")&df2$Enh.class%in%c("ME_nonprim","Global"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("HUES8_PKO_TET2_rescue")&df2$Enh.class%in%c("ME_prim","Global"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("HUES8_PKO_TET2_rescue")&df2$Enh.class%in%c("ME_nonprim","Global"),]))


summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("HUES8_DKO_P6")&df2$Enh.class%in%c("NPC_nonprim","NPC_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("HUES8_PKO_TET1s_rescue")&df2$Enh.class%in%c("NPC_nonprim","NPC_prim"),]))



####Stat test for human somatic tissue methyl####
enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(enh.lookup)=c("ID","Enh.class")
df=read.delim("./Data/human/DNA_methyl/Roadmap_2013_WGBS_v2.txt")[,c(1,13:30)]
colnames(df)[1]="ID"
df=merge(df,enh.lookup,by="ID")
df=pivot_longer(df[,c(20,2:19)],cols=-1,names_to = "Tissue",values_to="Methyl")
global=read.delim("./Data/human/DNA_methyl/Roadmap_2013_WGBS_Global.txt")[,c(7,13:30)]
global$ID=rep("Global")
colnames(global)[1]="Enh.class"
global=pivot_longer(global,cols=-1,names_to = "Tissue",values_to="Methyl")
df2=rbind(df,global)
enh.plot=c("Global","ME_nonprim","ME_prim","NPC_nonprim","NPC_prim","7W_brain_nonprim","7W_brain_prim")
tissue.plot=c("Adrenal","Aorta","Bladder","Fat","Liver","Lung","Thymus")
df2$Tissue=factor(df2$Tissue,levels=tissue.plot)
df2$Enh.class=factor(df2$Enh.class,levels=enh.plot)

summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Adrenal")&df2$Enh.class%in%c("ME_nonprim","ME_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Aorta")&df2$Enh.class%in%c("ME_nonprim","ME_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Bladder")&df2$Enh.class%in%c("ME_nonprim","ME_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Fat")&df2$Enh.class%in%c("ME_nonprim","ME_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Liver")&df2$Enh.class%in%c("ME_nonprim","ME_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Lung")&df2$Enh.class%in%c("ME_nonprim","ME_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Thymus")&df2$Enh.class%in%c("ME_nonprim","ME_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Enh.class%in%c("ME_nonprim","ME_prim"),]))


summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Adrenal")&df2$Enh.class%in%c("NPC_nonprim","NPC_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Aorta")&df2$Enh.class%in%c("NPC_nonprim","NPC_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Bladder")&df2$Enh.class%in%c("NPC_nonprim","NPC_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Fat")&df2$Enh.class%in%c("NPC_nonprim","NPC_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Liver")&df2$Enh.class%in%c("NPC_nonprim","NPC_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Lung")&df2$Enh.class%in%c("NPC_nonprim","NPC_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Thymus")&df2$Enh.class%in%c("NPC_nonprim","NPC_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Enh.class%in%c("NPC_nonprim","NPC_prim"),]))

summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Adrenal")&df2$Enh.class%in%c("7W_brain_nonprim","7W_brain_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Aorta")&df2$Enh.class%in%c("7W_brain_nonprim","7W_brain_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Bladder")&df2$Enh.class%in%c("7W_brain_nonprim","7W_brain_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Fat")&df2$Enh.class%in%c("7W_brain_nonprim","7W_brain_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Liver")&df2$Enh.class%in%c("7W_brain_nonprim","7W_brain_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Lung")&df2$Enh.class%in%c("7W_brain_nonprim","7W_brain_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Tissue%in%c("Thymus")&df2$Enh.class%in%c("7W_brain_nonprim","7W_brain_prim"),]))
summary(aov(Methyl~Enh.class,df2[df2$Enh.class%in%c("7W_brain_nonprim","7W_brain_prim"),]))



####Stat tests for gastrulation scNMT methylation plot####
enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(enh.lookup)=c("ID","Enh.class")
df=read.delim("./Data/human/DNA_methyl/Roadmap_2013_WGBS_v2.txt")[,c(1,13:30)]
head(df)
df2=read.delim("../single_cell/scnmt_gastrulation-develop/metacc/boxplots_enhancers/met_values2.txt")[,c(2:4)]
head(df2)
colnames(df2)=c("Tissue","Enh.class","Methyl")
df3=read.delim("../single_cell/scnmt_gastrulation-develop/metacc/boxplots_enhancers/Global_met_val.txt")[,c(4,5)]
colnames(df3)=c("Methyl","Tissue")
df3$Enh.class="Global"
df3=df3[,colnames(df2)]
df=rbind(df2,df3)
df$Tissue=str_replace(df$Tissue,"_"," ")

enh.plot=c("Global","mECT_prim","mECT_nonprim","mEND_prim","mEND_nonprim","mMES_prim","mMES_nonprim","forebrain_prim","forebrain_nonprim")
tissue.plot=c("E6.5 Epiblast","E7.5 Ectoderm","E7.5 Mesoderm","E7.5 Endoderm")
df=df[df$Tissue%in%tissue.plot,]
df=df[df$Enh.class%in%enh.plot,]
df$Tissue=factor(df$Tissue,levels=tissue.plot)
df$Enh.class=factor(df$Enh.class,levels=enh.plot)

summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E6.5 Epiblast")&df$Enh.class%in%c("Global","mECT_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E6.5 Epiblast")&df$Enh.class%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E6.5 Epiblast")&df$Enh.class%in%c("mEND_prim","mEND_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E6.5 Epiblast")&df$Enh.class%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E6.5 Epiblast")&df$Enh.class%in%c("forebrain_prim","forebrain_nonprim"),]))

summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Ectoderm")&df$Enh.class%in%c("Global","mECT_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Ectoderm")&df$Enh.class%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Ectoderm")&df$Enh.class%in%c("mEND_prim","mEND_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Ectoderm")&df$Enh.class%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Ectoderm")&df$Enh.class%in%c("forebrain_prim","forebrain_nonprim"),]))

summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Mesoderm")&df$Enh.class%in%c("Global","mECT_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Mesoderm")&df$Enh.class%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Mesoderm")&df$Enh.class%in%c("mEND_prim","mEND_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Mesoderm")&df$Enh.class%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Mesoderm")&df$Enh.class%in%c("forebrain_prim","forebrain_nonprim"),]))

summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Endoderm")&df$Enh.class%in%c("Global","mECT_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Endoderm")&df$Enh.class%in%c("mECT_prim","mECT_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Endoderm")&df$Enh.class%in%c("mEND_prim","mEND_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Endoderm")&df$Enh.class%in%c("mMES_prim","mMES_nonprim"),]))
summary(aov(Methyl~Enh.class,df[df$Tissue%in%c("E7.5 Endoderm")&df$Enh.class%in%c("forebrain_prim","forebrain_nonprim"),]))
