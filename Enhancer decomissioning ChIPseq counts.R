#is h3k4me1 lost upon differentiation into alternative lineages

setwd("G:/Headstone/Enh_priming/")

poised.lookup=read.delim("./Rebuttal/Annotations/Human_Specific_Poised_enhancers.txt",h=F)[,4:5]
enh.lookup=read.delim("./Annotations/Human/Comb_Enh_class_lookup_v2.txt",h=F)
colnames(enh.lookup)=colnames(poised.lookup)=c("Probe","Enh.class")
head(enh.lookup)
enh.lookup=rbind(enh.lookup,poised.lookup)
summary(as.factor(enh.lookup$Enh.class))

rpkm=read.delim("./Rebuttal/Data/Xie_et_al_ChIP_RPKM_500bp_core.txt")
rpkm=rpkm[,c(1,13:24)]

m.rpkm=merge(rpkm,enh.lookup,by="Probe")
library(dplyr)
library(tidyverse)
library(ggplot2)
help(pivot_longer)
head(m.rpkm)
df=m.rpkm %>% pivot_longer(cols=H1_H3K27ac:ME_H3K27me3, names_to="Tissue",values_to = "RPKM")
head(df)
unique(df$Enh.class)
to.plot=c("brain_poised","7W_brain_prim","7W_brain_nonprim",
          "liver_poised","Fetal_liver_prim","Fetal_liver_nonprim",
          "kidney_poised","Fetal_kidney_prim","Fetal_kidney_nonprim",
          "NPC_poised","NPC_prim","NPC_nonprim",
          "ME_poised","ME_prim","ME_nonprim",
          "MSC_poised","MSC_prim","MSC_nonprim")
{
#ggplot(df[df$Enh.class%in%to.plot,],aes(x=Tissue,y=RPKM,fill=Enh.class))+geom_boxplot(outlier.shape = NA)
#ggplot(df[df$Enh.class%in%to.plot&df$Tissue=="NPC_H3K4me1",],aes(x=Tissue,y=RPKM,fill=Enh.class))+geom_boxplot(outlier.shape = NA)

#df.sub=df[df$Enh.class%in%to.plot,]
#df.sub$Enh.class=factor(df.sub$Enh.class,levels=to.plot)

#marks=unique(df.sub$Tissue)
#for(i in marks){
#plot(ggplot(df.sub[df.sub$Tissue==i,],aes(x=Enh.class,y=log(RPKM+0.01),fill=Enh.class))+geom_boxplot(outlier.shape = NA)+ggtitle(i))
#}
}
global=read.delim("./Rebuttal/Data/Xie_et_al_ChIP_RPKM_2M_random_ctrl.txt")
head(global)
exon.instances=global$Feature!="null"
global.out=global[!exon.instances,]
exon.counts=global[exon.instances,]

g.df=global.out[,c(1,13:24)]%>%pivot_longer(cols=H1_H3K27ac:MSC_H3K27me3, names_to="Tissue",values_to = "RPKM")
g.df$Enh.class=rep("Global")
e.df=exon.counts[,c(1,13:24)]%>%pivot_longer(cols=H1_H3K27ac:MSC_H3K27me3, names_to="Tissue",values_to = "RPKM")
e.df$Enh.class=rep("Exon")

g.df=g.df[,colnames(df)]
e.df=e.df[,colnames(df)]
tog.df=rbind(df,g.df)
tog.df=rbind(tog.df,e.df)

to.plot2=c("Global","Exon",
          "brain_poised","7W_brain_prim","7W_brain_nonprim",
          "liver_poised","Fetal_liver_prim","Fetal_liver_nonprim",
          "kidney_poised","Fetal_kidney_prim","Fetal_kidney_nonprim",
          "NPC_poised","NPC_prim","NPC_nonprim",
          "ME_poised","ME_prim","ME_nonprim",
          "MSC_poised","MSC_prim","MSC_nonprim")
tog.df=tog.df[tog.df$Enh.class%in%to.plot2,]
tog.df$Enh.class=factor(tog.df$Enh.class,levels=to.plot2)
marks=unique(tog.df$Tissue)
#for(i in marks){
#  plot(ggplot(tog.df[tog.df$Tissue==i,],aes(x=Enh.class,y=log(RPKM+0.01),fill=Enh.class))+geom_boxplot(outlier.shape = NA)+ggtitle(i))
#}
#for(i in marks){
#  plot(ggplot(tog.df[tog.df$Tissue==i,],aes(x=Enh.class,y=RPKM,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+ggtitle(i))
#}

#i="ME_H3K27me3"
#plot(ggplot(tog.df[tog.df$Tissue==i,],aes(x=Enh.class,y=RPKM,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+ggtitle(i))
#plot(ggplot(tog.df[tog.df$Tissue==i,],aes(x=Enh.class,y=RPKM,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+ggtitle(i))+coord_cartesian(ylim=c(0,18))

uppers=c(10,10,30,
         10,10,10,
         10,10,10,
         10,10,20)
pdf("./Rebuttal/Plots/ChIP_seq_RPKM_plots.pdf")
for(i in 1:length(marks)){
  plot(ggplot(tog.df[tog.df$Tissue==marks[i],],aes(x=Enh.class,y=RPKM,fill=Enh.class))+geom_boxplot(outlier.shape = NA)+ggtitle(marks[i])+coord_cartesian(ylim=c(0,uppers[i]))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
}
dev.off()
