setwd("G:/Headstone/Enh_priming/")
library(tidyverse)
library(networkD3)
library(GenomicRanges)
library(rtracklayer)
library(htmlwidgets)
library(webshot)



enh.class="MSC"
class.spec.col=(5)
#1-Probe:2-Other:3-H1:4-ME:5-MSC:6-NPC
{

xie.enh=read.delim("./Annotations/Human/Xie_et_al_2013_coord_hg38.bed",h=F)
class.all=nrow(unique(xie.enh[grep(paste0(enh.class,"_"),xie.enh$V4),1:3]))

tss.removed.xie=read.delim("./Annotations/Human/Xie_enh_coord_hg38.txt",h=F)
class.tss=unique(tss.removed.xie[grep(paste0(enh.class,"_"),tss.removed.xie$V5),4])


#loading the RPKM values
xie.rpkms=read.delim("./Data/human/ChIP/Xie_Enh_1_5kb_ChIP_logRPKM_v2.txt")[,c(1,13:18)]
fetal.rpkms=read.delim("./Data/human/ChIP/Fetal_Enh_1_5kb_ChIP_logRPKM.txt")[,c(1,13:18)]
brain.rpkms=read.delim("./Data/human/ChIP/7W_Brain_Enh_1_5kb_ChIP_logRPKM_v2.txt")[,c(1,13:18)]
tog.rpkms=do.call(rbind,list(xie.rpkms,fetal.rpkms,brain.rpkms))

set.k4me1=1.2
set.k27ac=2.97
set.k27me3=2

all.prim=tog.rpkms$Probe[tog.rpkms$H1_K27ac<set.k27ac&tog.rpkms$H1_K4me1>set.k4me1&tog.rpkms$H1_K27me3<set.k27me3]
all.poised=tog.rpkms$Probe[tog.rpkms$H1_K27ac<set.k27ac&tog.rpkms$H1_K27me3>=set.k27me3]

class.prim=all.prim[all.prim%in%class.tss]
class.poised=all.poised[all.poised%in%class.tss]



tss.removed.xie$V1=gsub("chr","",tss.removed.xie$V1)
peaks.h1=import.bed("./Annotations/Human/MACS/H1_H3K27ac.bed")
peaks.ME=import.bed("./Annotations/Human/MACS/ME_H3K27ac.bed")
peaks.MSC=import.bed("./Annotations/Human/MACS/MSC_H3K27ac.bed")
peaks.NPC=import.bed("./Annotations/Human/MACS/NPC_H3K27ac.bed")

xie.class=enh.class
{
  xie.sub=unique(tss.removed.xie[grep(pattern = paste0(xie.class,"_enriched"),tss.removed.xie$V5),1:4])
  class.enriched=unique(xie.sub$V4)
  xie.sub=unique(tss.removed.xie[grep(pattern = paste0(xie.class),tss.removed.xie$V5),1:4])
  xie.other=unique(tss.removed.xie[!grepl(pattern = xie.class,tss.removed.xie$V5),1:4])
  tmp.sub=tempfile()
  tmp.other=tempfile()
  write.table(xie.sub,tmp.sub,sep="\t",col.names = F,row.names = F,quote = F)
  write.table(xie.other,tmp.other,sep="\t",col.names = F,row.names = F,quote = F)
  sub.coord=import.bed(tmp.sub)
  other.coord=import.bed(tmp.other)
  file.remove(tmp.sub)
  file.remove(tmp.other)
  rm(xie.sub,xie.other)
  sub.coord=sub.coord%>%GenomicRanges::resize(width=100,fix="center")
  other.coord=other.coord%>%GenomicRanges::resize(width=100,fix="center")
}
#sub.coord=sub.coord[sub.coord %outside% other.coord]
fo.other=findOverlaps(sub.coord,other.coord)
sub.coord=sub.coord%>%GenomicRanges::resize(width=1500,fix="center")
fo.h1=findOverlaps(sub.coord,peaks.h1)
fo.ME=findOverlaps(sub.coord,peaks.ME)
fo.MSC=findOverlaps(sub.coord,peaks.MSC)
fo.NPC=findOverlaps(sub.coord,peaks.NPC)

df=data.frame("ID"=sub.coord$name,
              "other"=rep(F),
              "H1"=rep(F),
              "ME"=rep(F),
              "MSC"=rep(F),
              "NPC"=rep(F))
df$other[queryHits(fo.other)]=T
df$H1[queryHits(fo.h1)]=T
df$ME[queryHits(fo.ME)]=T
df$MSC[queryHits(fo.MSC)]=T
df$NPC[queryHits(fo.NPC)]=T


class.cols=c(2:6)
class.cols=class.cols[!class.cols%in%class.spec.col]


#1-Probe:2-Other:3-H1:4-ME:5-MSC:6-NPC

embryonic.col=3
other.cols= (4:6)
other.cols=other.cols[!other.cols%in%class.spec.col]

class.spec=df$ID[df[,embryonic.col]==F&rowSums(df[,other.cols])==0]
hESC=df$ID[df[,embryonic.col]==T&rowSums(df[,other.cols])==0]
Other.lineage=df$ID[df[,embryonic.col]==F&rowSums(df[,other.cols])>0]
Emrbyonic.and.other=df$ID[df[,embryonic.col]==T&rowSums(df[,other.cols])>0]


class.other=df$ID[df$other]
class.noother=df$ID[df$other==F]

class.prim.enr=class.prim[class.prim%in%class.enriched]
class.prim.nonenr=class.prim[!class.prim%in%class.enriched]
class.prim.enr.noother=class.prim.enr[!class.prim.enr%in%df$ID[df$other]]
class.prim.enr.other=class.prim.enr[class.prim.enr%in%df$ID[df$other]]


class.poised.enr=class.poised[class.poised%in%class.enriched]
class.poised.nonenr=class.poised[!class.poised%in%class.enriched]
class.poised.enr.noother=class.poised.enr[!class.poised.enr%in%df$ID[df$other]]
class.poised.enr.other=class.poised.enr[class.poised.enr%in%df$ID[df$other]]


#Sankey functions
#initial primed poised
{
  data <- matrix(rep(0), 6, 6)
  colnames(data) = rownames(data) = c("All Enhancers","Enhancers of class","Not overlap TSS","Overlap TSS", "Primed", "Poised")
  data[2,3]=length(class.tss)
  data[2,4]=class.all-length(class.tss)
  data[3,5]=length(class.tss[class.tss%in%all.prim])
  data[3,6]=length(class.tss[class.tss%in%all.poised])
  links <- data %>% 
    as.data.frame() %>% 
    rownames_to_column(var="source") %>% 
    gather(key="target", value="value", -1) %>%
    filter(value != 0)
  
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  links
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     sinksRight=FALSE)
  p
  saveWidget(p,file="./Rebuttal/sankey_test.html")
  webshot("./Rebuttal/sankey_test.html" , paste0("./Sankey_plots/",enh.class,"_primed_vs_poised.pdf"), delay = 0.2)
  write.table(links,paste0("./Sankey_plots/",enh.class,"_primed_vs_poised.txt"),col.names = T,row.names = F,quote = F)
  
}
#primed specific non specific
{
  
  data <- matrix(rep(0), 9, 9)
  colnames(data) = rownames(data) = c("Primed","Enriched","Not Enriched","No Overlap" ,"Overlap Other", "No Overlap2", "Overlap Other2","Primed Specific","Primed Non-specific")
  
  data[1,2]=length(class.prim.enr)
  data[1,3]=length(class.prim.nonenr)
  data[2,4]=length(class.prim.enr.noother)
  data[2,5]=length(class.prim.enr.other)
  data[3,6]=length(class.prim.nonenr[class.prim.nonenr%in%class.noother])
  data[3,7]=length(class.prim.nonenr[class.prim.nonenr%in%class.other])
  data[4,8]=length(class.prim.enr.noother)
  data[5,9]=data[2,5]
  data[6,9]=data[3,6]
  data[7,9]=data[3,7]
  links <- data %>% 
    as.data.frame() %>% 
    rownames_to_column(var="source") %>% 
    gather(key="target", value="value", -1) %>%
    filter(value != 0)
  
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  links
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     sinksRight=FALSE)
  p
  saveWidget(p,file="./Rebuttal/sankey_test.html")
  webshot("./Rebuttal/sankey_test.html" , paste0("./Sankey_plots/",enh.class,"_primed_specific.pdf"), delay = 0.2)
  write.table(links,paste0("./Sankey_plots/",enh.class,"_primed_specific.txt"),col.names = T,row.names = F,quote = F)
  primed.specific=class.prim.enr.noother
  primed.nonspecific=c(class.prim.enr.other,class.prim.nonenr)
}
#primed additional filtering
{
  
  data <- matrix(rep(0), 10, 10)
  colnames(data) = rownames(data) = c("Specific","Non-Specific","Tissue only","+Embryonic" ,"+other", "+Embryonic+other","Tissue only2","+Embryonic2" ,"+other2", "+Embryonic+other2")
  
  data[1,3]=length(primed.specific[primed.specific%in%class.spec])
  data[1,4]=length(primed.specific[primed.specific%in%hESC])
  data[1,5]=length(primed.specific[primed.specific%in%Other.lineage])
  data[1,6]=length(primed.specific[primed.specific%in%Emrbyonic.and.other])
  
  data[2,7]=length(primed.nonspecific[primed.nonspecific%in%class.spec])
  data[2,8]=length(primed.nonspecific[primed.nonspecific%in%hESC])
  data[2,9]=length(primed.nonspecific[primed.nonspecific%in%Other.lineage])
  data[2,10]=length(primed.nonspecific[primed.nonspecific%in%Emrbyonic.and.other])
  

  links <- data %>% 
    as.data.frame() %>% 
    rownames_to_column(var="source") %>% 
    gather(key="target", value="value", -1) %>%
    filter(value != 0)
  
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  links
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     sinksRight=FALSE)
  p
  saveWidget(p,file="./Rebuttal/sankey_test.html")
  webshot("./Rebuttal/sankey_test.html" , paste0("./Sankey_plots/",enh.class,"_primed_additional.pdf"), delay = 0.2)
  write.table(links,paste0("./Sankey_plots/",enh.class,"_primed_additional.txt"),col.names = T,row.names = F,quote = F)
  
}
#poised specific non specific
{
  
  data <- matrix(rep(0), 9, 9)
  colnames(data) = rownames(data) = c("Poised","Enriched","Not Enriched","No Overlap" ,"Overlap Other", "No Overlap2", "Overlap Other2","Poised Specific","Poised Non-specific")
  
  data[1,2]=length(class.poised.enr)
  data[1,3]=length(class.poised.nonenr)
  data[2,4]=length(class.poised.enr.noother)
  data[2,5]=length(class.poised.enr.other)
  data[3,6]=length(class.poised.nonenr[class.poised.nonenr%in%class.noother])
  data[3,7]=length(class.poised.nonenr[class.poised.nonenr%in%class.other])
  data[4,8]=length(class.poised.enr.noother)
  data[5,9]=data[2,5]
  data[6,9]=data[3,6]
  data[7,9]=data[3,7]
  links <- data %>% 
    as.data.frame() %>% 
    rownames_to_column(var="source") %>% 
    gather(key="target", value="value", -1) %>%
    filter(value != 0)
  
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  links
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     sinksRight=FALSE)
  p
  poised.specific=class.poised.enr.noother
  poised.nonspecific=c(class.poised.enr.other,class.poised.nonenr)
  saveWidget(p,file="./Rebuttal/sankey_test.html")
  webshot("./Rebuttal/sankey_test.html" , paste0("./Sankey_plots/",enh.class,"_poised_specific.pdf"), delay = 0.2)
  write.table(links,paste0("./Sankey_plots/",enh.class,"_poised_specific.txt"),col.names = T,row.names = F,quote = F)
  
}
#poised additional filtering
{
  
  data <- matrix(rep(0), 10, 10)
  colnames(data) = rownames(data) = c("Specific","Non-Specific","Tissue only","+Embryonic" ,"+other", "+Embryonic+other","Tissue only2","+Embryonic2" ,"+other2", "+Embryonic+other2")
  
  data[1,3]=length(poised.specific[poised.specific%in%class.spec])
  data[1,4]=length(poised.specific[poised.specific%in%hESC])
  data[1,5]=length(poised.specific[poised.specific%in%Other.lineage])
  data[1,6]=length(poised.specific[poised.specific%in%Emrbyonic.and.other])
  
  data[2,7]=length(poised.nonspecific[poised.nonspecific%in%class.spec])
  data[2,8]=length(poised.nonspecific[poised.nonspecific%in%hESC])
  data[2,9]=length(poised.nonspecific[poised.nonspecific%in%Other.lineage])
  data[2,10]=length(poised.nonspecific[poised.nonspecific%in%Emrbyonic.and.other])
  
  
  links <- data %>% 
    as.data.frame() %>% 
    rownames_to_column(var="source") %>% 
    gather(key="target", value="value", -1) %>%
    filter(value != 0)
  
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  links
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     sinksRight=FALSE)
  p
  saveWidget(p,file="./Rebuttal/sankey_test.html")
  webshot("./Rebuttal/sankey_test.html" , paste0("./Sankey_plots/",enh.class,"_poised_additional.pdf"), delay = 0.2)
  write.table(links,paste0("./Sankey_plots/",enh.class,"_poised_additional.txt"),col.names = T,row.names = F,quote = F)

}
#primed full
{
  data <- matrix(rep(0), 13, 13)
  colnames(data) = rownames(data) = c("Enhancers of class","Not overlap TSS","Overlap TSS", "Primed", "Poised","Enriched","Not Enriched","No Overlap" ,"Overlap Other","Tissue only","+Embryonic" ,"+other", "+Embryonic+other")
  
  data[1,2]=length(class.tss)
  data[1,3]=class.all-length(class.tss)
  
  data[2,4]=length(class.tss[class.tss%in%all.prim])
  data[2,5]=length(class.tss[class.tss%in%all.poised])
  
  data[4,6]=length(class.prim.enr)
  data[4,7]=length(class.prim.nonenr)
  
  data[6,8]=length(class.prim.enr.noother)
  data[6,9]=length(class.prim.enr.other)
  
  data[8,10]=length(primed.specific[primed.specific%in%class.spec])
  data[8,11]=length(primed.specific[primed.specific%in%hESC])
  data[8,12]=length(primed.specific[primed.specific%in%Other.lineage])
  data[8,13]=length(primed.specific[primed.specific%in%Emrbyonic.and.other])
  
  
  links <- data %>% 
    as.data.frame() %>% 
    rownames_to_column(var="source") %>% 
    gather(key="target", value="value", -1) %>%
    filter(value != 0)
  
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  links
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     sinksRight=FALSE)
  p
  saveWidget(p,file="./Rebuttal/sankey_test.html")
  webshot("./Rebuttal/sankey_test.html" , paste0("./Sankey_plots/",enh.class,"_primed_full.pdf"), delay = 0.2)
  write.table(links,paste0("./Sankey_plots/",enh.class,"_primed_full.txt"),col.names = T,row.names = F,quote = F)
  
  
}
#poised full
{
  data <- matrix(rep(0), 13, 13)
  colnames(data) = rownames(data) = c("Enhancers of class","Not overlap TSS","Overlap TSS", "Primed", "Poised","Enriched","Not Enriched","No Overlap" ,"Overlap Other","Tissue only","+Embryonic" ,"+other", "+Embryonic+other")
  
  data[1,2]=length(class.tss)
  data[1,3]=class.all-length(class.tss)
  
  data[2,4]=length(class.tss[class.tss%in%all.prim])
  data[2,5]=length(class.tss[class.tss%in%all.poised])
  
  data[5,6]=length(class.poised.enr)
  data[5,7]=length(class.poised.nonenr)
  
  data[6,8]=length(class.poised.enr.noother)
  data[6,9]=length(class.poised.enr.other)
  
  data[8,10]=length(poised.specific[poised.specific%in%class.spec])
  data[8,11]=length(poised.specific[poised.specific%in%hESC])
  data[8,12]=length(poised.specific[poised.specific%in%Other.lineage])
  data[8,13]=length(poised.specific[poised.specific%in%Emrbyonic.and.other])
  
  
  links <- data %>% 
    as.data.frame() %>% 
    rownames_to_column(var="source") %>% 
    gather(key="target", value="value", -1) %>%
    filter(value != 0)
  
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  links
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     sinksRight=FALSE)
  p

  saveWidget(p,file="./Rebuttal/sankey_test.html")
  webshot("./Rebuttal/sankey_test.html" , paste0("./Sankey_plots/",enh.class,"_poised_full.pdf"), delay = 0.2)
  write.table(links,paste0("./Sankey_plots/",enh.class,"_poised_full.txt"),col.names = T,row.names = F,quote = F)
  
}
}

####Brain####
enh.class="7W_brain"
import.remove.chr=function(file){
  x=read.delim(file,h=F)
  x[,1]=gsub("chr","",x[,1])
  tmp=tempfile()
  write.table(x,tmp,sep="\t",col.names = F,row.names = F,quote = F)
  gr=import.bed(tmp)
  file.remove(tmp)
  return(gr)
}
TSS=import.remove.chr("./Rebuttal/Human_ENSEMBL_TSS.bed")
TSS=TSS%>%GenomicRanges::resize(width=200,fix="center")

##7 week fetal brain peaks
brain.peaks=import.remove.chr("./Annotations/Human/MACS/7pcw_enh_hg38.bed")
num.brain.peaks=length(width(brain.peaks))

brain.non.prom=brain.peaks[brain.peaks%outside%TSS]
num.brain.tss=length(width(brain.non.prom))


brain.non.prom=import.remove.chr("./Annotations/Human/Enhancers/7W_Brain_Enh.txt")
brain.ids=brain.non.prom$name

xie.rpkms=read.delim("./Data/human/ChIP/Xie_Enh_1_5kb_ChIP_logRPKM_v2.txt")[,c(1,13:18)]
fetal.rpkms=read.delim("./Data/human/ChIP/Fetal_Enh_1_5kb_ChIP_logRPKM.txt")[,c(1,13:18)]
brain.rpkms=read.delim("./Data/human/ChIP/7W_Brain_Enh_1_5kb_ChIP_logRPKM_v2.txt")[,c(1,13:18)]
tog.rpkms=do.call(rbind,list(xie.rpkms,fetal.rpkms,brain.rpkms))


set.k4me1=1.2
set.k27ac=2.97
set.k27me3=2

all.prim=tog.rpkms$Probe[tog.rpkms$H1_K27ac<set.k27ac&tog.rpkms$H1_K4me1>set.k4me1&tog.rpkms$H1_K27me3<set.k27me3]
all.poised=tog.rpkms$Probe[tog.rpkms$H1_K27ac<set.k27ac&tog.rpkms$H1_K27me3>=set.k27me3]


brain.prim=brain.ids[brain.ids%in%all.prim]
brain.poised=brain.ids[brain.ids%in%all.poised]

tss.removed.xie=unique(read.delim("./Annotations/Human/Xie_enh_coord_hg38.txt",h=F)[,1:4])
tss.removed.xie$V1=gsub("chr","",tss.removed.xie$V1)
tmp=tempfile()
write.table(tss.removed.xie,tmp,sep="\t",col.names = F,row.names = F,quote = F)
xie.coord=import.bed(tmp)
file.remove(tmp)
rm(tss.removed.xie)

peaks.h1=import.bed("./Annotations/Human/MACS/H1_H3K27ac.bed")
peaks.ME=import.bed("./Annotations/Human/MACS/ME_H3K27ac.bed")
peaks.MSC=import.bed("./Annotations/Human/MACS/MSC_H3K27ac.bed")
peaks.NPC=import.bed("./Annotations/Human/MACS/NPC_H3K27ac.bed")

##Getting lineage specific 7pcw Brain enh


fo.other=findOverlaps(brain.non.prom,xie.coord)
brain.non.prom=brain.non.prom%>%GenomicRanges::resize(width=1500,fix="center")
fo.h1=findOverlaps(brain.non.prom,peaks.h1)
fo.ME=findOverlaps(brain.non.prom,peaks.ME)
fo.MSC=findOverlaps(brain.non.prom,peaks.MSC)
fo.NPC=findOverlaps(brain.non.prom,peaks.NPC)

df=data.frame("ID"=brain.non.prom$name,
              "other"=rep(F),
              "H1"=rep(F),
              "ME"=rep(F),
              "MSC"=rep(F),
              "NPC"=rep(F))
df$other[queryHits(fo.other)]=T
df$H1[queryHits(fo.h1)]=T
df$ME[queryHits(fo.ME)]=T
df$MSC[queryHits(fo.MSC)]=T
df$NPC[queryHits(fo.NPC)]=T

embryonic.col=3
other.cols= (4:6)

overlap.xie=df$ID[df$other]

class.spec=df$ID[df[,embryonic.col]==F&rowSums(df[,other.cols])==0]
hESC=df$ID[df[,embryonic.col]==T&rowSums(df[,other.cols])==0]
Other.lineage=df$ID[df[,embryonic.col]==F&rowSums(df[,other.cols])>0]
Emrbyonic.and.other=df$ID[df[,embryonic.col]==T&rowSums(df[,other.cols])>0]

#primed full
{
  data <- matrix(rep(0), 11, 11)
  colnames(data) = rownames(data) = c("Number of Brain H3K27ac peaks","Not overlap TSS","Overlap TSS", "Primed", "Poised","No Overlap" ,"Overlap Early Dev Enh","Tissue only","+Embryonic" ,"+Early Dev", "+Embryonic+Early Dev")
  
  data[1,2]=num.brain.tss
  data[1,3]=num.brain.peaks-num.brain.tss
  
  data[2,4]=length(brain.ids[brain.ids%in%all.prim])
  data[2,5]=length(brain.ids[brain.ids%in%all.poised])
  
  data[4,6]=length(brain.prim[!brain.prim%in%overlap.xie])
  data[4,7]=length(brain.prim[brain.prim%in%overlap.xie])
  primed.specific=brain.prim[!brain.prim%in%overlap.xie]
  data[6,8]=length(primed.specific[primed.specific%in%class.spec])
  data[6,9]=length(primed.specific[primed.specific%in%hESC])
  data[6,10]=length(primed.specific[primed.specific%in%Other.lineage])
  data[6,11]=length(primed.specific[primed.specific%in%Emrbyonic.and.other])
  
  
  
  links <- data %>% 
    as.data.frame() %>% 
    rownames_to_column(var="source") %>% 
    gather(key="target", value="value", -1) %>%
    filter(value != 0)
  
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  links
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     sinksRight=FALSE)
  p
  saveWidget(p,file="./Rebuttal/sankey_test.html")
  webshot("./Rebuttal/sankey_test.html" , paste0("./Sankey_plots/",enh.class,"_primed_full.pdf"), delay = 0.2)
  write.table(links,paste0("./Sankey_plots/",enh.class,"_primed_full.txt"),col.names = T,row.names = F,quote = F)
  
  
}
#poised full
{
  data <- matrix(rep(0), 11, 11)
  colnames(data) = rownames(data) = c("Number of Brain H3K27ac peaks","Not overlap TSS","Overlap TSS", "Primed", "Poised","No Overlap" ,"Overlap Early Dev Enh","Tissue only","+Embryonic" ,"+Early Dev", "+Embryonic+Early Dev")
  
  data[1,2]=num.brain.tss
  data[1,3]=num.brain.peaks-num.brain.tss
  
  data[2,4]=length(brain.ids[brain.ids%in%all.prim])
  data[2,5]=length(brain.ids[brain.ids%in%all.poised])
  
  data[5,6]=length(brain.poised[!brain.poised%in%overlap.xie])
  data[5,7]=length(brain.poised[brain.poised%in%overlap.xie])
  poised.specific=brain.poised[!brain.poised%in%overlap.xie]
  data[6,8]=length(poised.specific[poised.specific%in%class.spec])
  data[6,9]=length(poised.specific[poised.specific%in%hESC])
  data[6,10]=length(poised.specific[poised.specific%in%Other.lineage])
  data[6,11]=length(poised.specific[poised.specific%in%Emrbyonic.and.other])
  
  
  links <- data %>% 
    as.data.frame() %>% 
    rownames_to_column(var="source") %>% 
    gather(key="target", value="value", -1) %>%
    filter(value != 0)
  
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  links
  # Make the Network
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     sinksRight=FALSE)
  p
  
  saveWidget(p,file="./Rebuttal/sankey_test.html")
  webshot("./Rebuttal/sankey_test.html" , paste0("./Sankey_plots/",enh.class,"_poised_full.pdf"), delay = 0.2)
  write.table(links,paste0("./Sankey_plots/",enh.class,"_poised_full.txt"),col.names = T,row.names = F,quote = F)
  
}
