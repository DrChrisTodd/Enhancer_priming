###Generating the enhancer subsets which will be used to plot the alluvials

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


plot.alluvial=function(enh.type){

all.ids=human.enh$ID[human.enh$Enh.class==enh.type]
rpkm.sub=tog.rpkms[tog.rpkms$Probe%in%all.ids,]
cutoff.non.prim=max(rpkm.sub$H1_K4me1[rpkm.sub$Probe%in%read.delim(paste0("./Annotations/Human/Enhancers/",annotation.files[grepl(paste0(enh.type,"_nonprim"),annotation.files)]),h=F)[,4]])
prim=rpkm.sub$Probe[rpkm.sub$H1_K27ac<set.k27ac&rpkm.sub$H1_K27me3<set.k27me3&rpkm.sub$H1_K4me1>set.k4me1]
poised=rpkm.sub$Probe[rpkm.sub$H1_K27ac<set.k27ac&rpkm.sub$H1_K27me3>=set.k27me3]
non.prim=rpkm.sub$Probe[rpkm.sub$H1_K27ac<set.k27ac&rpkm.sub$H1_K27me3<set.k27me3&rpkm.sub$H1_K4me1<=cutoff.non.prim]
active=rpkm.sub$Probe[rpkm.sub$H1_K27ac>=set.k27ac&rpkm.sub$H1_K27me3<set.k27me3]
intermediate=all.ids[!all.ids%in%c(prim,poised,non.prim,active)]

specific.ids=unique(comb.lookup$ID[grepl(enh.type,comb.lookup$Enh.class)])
non.specific.ids=all.ids[!all.ids%in%specific.ids]

{
df.rows=list()
df.rows[[1]]=data.frame(Lineage="Specific",class="active",feq=sum(specific.ids%in%active))
df.rows[[2]]=data.frame(Lineage="Specific",class="poised",feq=sum(specific.ids%in%poised))
df.rows[[3]]=data.frame(Lineage="Specific",class="primed",feq=sum(specific.ids%in%prim))
df.rows[[4]]=data.frame(Lineage="Specific",class="intermediate",feq=sum(specific.ids%in%intermediate))
df.rows[[5]]=data.frame(Lineage="Specific",class="nonprimed",feq=sum(specific.ids%in%non.prim))
df.rows[[6]]=data.frame(Lineage="Nonspecific",class="active",feq=sum(non.specific.ids%in%active))
df.rows[[7]]=data.frame(Lineage="Nonspecific",class="poised",feq=sum(non.specific.ids%in%poised))
df.rows[[8]]=data.frame(Lineage="Nonspecific",class="primed",feq=sum(non.specific.ids%in%prim))
df.rows[[9]]=data.frame(Lineage="Nonspecific",class="intermediate",feq=sum(non.specific.ids%in%intermediate))
df.rows[[10]]=data.frame(Lineage="Nonspecific",class="nonprimed",feq=sum(non.specific.ids%in%non.prim))
tog.df=do.call(rbind,df.rows)
tog.df=tog.df[,c(2,1,3)]
}
tog.df$class=factor(tog.df$class,levels=c("active","poised","primed","intermediate","nonprimed"))
tog.df$Lineage=factor(tog.df$Lineage,levels=c("Specific","Nonspecific"))

is_alluvia_form(tog.df,axes=1:3)

gg=ggplot(tog.df,aes(y=feq,axis1=class,axis2=Lineage))+
  geom_alluvium(aes(fill=class),width=1/24)+
  geom_stratum(width=1/12,fill="black",color="grey")+
  geom_label(stat="stratum",aes(label=after_stat(stratum)))+
  scale_x_discrete(limits=c("class","Lineage"),expand=c(.05,.05))+
  scale_fill_brewer(type="qual",palette="Set1")+
  ggtitle(paste0(enh.type," enhancer type breakdown"))+theme_bw()
plot(gg)
}

pdf("./Plots/Alluvial_enh_groups/Human_enh_group_alluvials.pdf")
plot.alluvial("NPC")
plot.alluvial("ME")
plot.alluvial("MSC")
plot.alluvial("brain")
plot.alluvial("liver")
plot.alluvial("kidney")
plot.alluvial("adrenal")
plot.alluvial("kidney")
plot.alluvial("ventricle")
plot.alluvial("lung")
plot.alluvial("stomach")
dev.off()


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

plot.alluvial=function(enh.type){
  
  all.ids=mouse.enh$ID[mouse.enh$Enh.class==enh.type]
  rpkm.sub=tog.rpkms[tog.rpkms$ID%in%all.ids,]
  cutoff.non.prim=max(rpkm.sub$K4me1[rpkm.sub$ID%in%comb.lookup$ID[comb.lookup$Enh.class==paste0(enh.type,"_nonprim")]])
  prim=rpkm.sub$ID[rpkm.sub$K27ac<set.k27ac&rpkm.sub$K27me3<set.k27me3&rpkm.sub$K4me1>set.k4me1]
  poised=rpkm.sub$ID[rpkm.sub$K27ac<set.k27ac&rpkm.sub$K27me3>=set.k27me3]
  non.prim=rpkm.sub$ID[rpkm.sub$K27ac<set.k27ac&rpkm.sub$K27me3<set.k27me3&rpkm.sub$K4me1<=cutoff.non.prim]
  active=rpkm.sub$ID[rpkm.sub$K27ac>=set.k27ac&rpkm.sub$K27me3<set.k27me3]
  intermediate=all.ids[!all.ids%in%c(prim,poised,non.prim,active)]
  
  specific.ids=unique(comb.lookup$ID[grepl(enh.type,comb.lookup$Enh.class)])
  non.specific.ids=all.ids[!all.ids%in%specific.ids]
  
  {
    df.rows=list()
    df.rows[[1]]=data.frame(Lineage="Specific",class="active",feq=sum(specific.ids%in%active))
    df.rows[[2]]=data.frame(Lineage="Specific",class="poised",feq=sum(specific.ids%in%poised))
    df.rows[[3]]=data.frame(Lineage="Specific",class="primed",feq=sum(specific.ids%in%prim))
    df.rows[[4]]=data.frame(Lineage="Specific",class="intermediate",feq=sum(specific.ids%in%intermediate))
    df.rows[[5]]=data.frame(Lineage="Specific",class="nonprimed",feq=sum(specific.ids%in%non.prim))
    df.rows[[6]]=data.frame(Lineage="Nonspecific",class="active",feq=sum(non.specific.ids%in%active))
    df.rows[[7]]=data.frame(Lineage="Nonspecific",class="poised",feq=sum(non.specific.ids%in%poised))
    df.rows[[8]]=data.frame(Lineage="Nonspecific",class="primed",feq=sum(non.specific.ids%in%prim))
    df.rows[[9]]=data.frame(Lineage="Nonspecific",class="intermediate",feq=sum(non.specific.ids%in%intermediate))
    df.rows[[10]]=data.frame(Lineage="Nonspecific",class="nonprimed",feq=sum(non.specific.ids%in%non.prim))
    tog.df=do.call(rbind,df.rows)
  }
  tog.df$class=factor(tog.df$class,levels=c("active","poised","primed","intermediate","nonprimed"))
  tog.df$Lineage=factor(tog.df$Lineage,levels=c("Specific","Nonspecific"))
  
  #is_alluvia_form(tog.df,axes=1:3)
  
  gg=ggplot(tog.df,aes(y=feq,axis1=class,axis2=Lineage))+
    geom_alluvium(aes(fill=class),width=1/24)+
    geom_stratum(width=1/12,fill="black",color="grey")+
    geom_label(stat="stratum",aes(label=after_stat(stratum)))+
    scale_x_discrete(limits=c("class","Lineage"),expand=c(.05,.05))+
    scale_fill_brewer(type="qual",palette="Set1")+
    ggtitle(paste0(enh.type," enhancer type breakdown"))+theme_bw()
  plot(gg)
}


pdf("./Plots/Alluvial_enh_groups/Mouse_enh_group_alluvials.pdf")
plot.alluvial("mECT")
plot.alluvial("mEND")
plot.alluvial("mMES")

plot.alluvial("forebrain")
plot.alluvial("heart")
plot.alluvial("hindbrain")
plot.alluvial("liver")
plot.alluvial("midbrain")
dev.off()

