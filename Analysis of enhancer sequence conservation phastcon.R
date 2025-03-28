library(data.table)
library(dplyr)
library(ggplot2)

species <- "mouse"
regions_human <- fread("~/subprojects/chris_paper/phastcons/Human_enhancers_1kb_reformat.hg38.phastCons20way.bed") %>% setnames(c("chr","start","end","name","class","conservation"))
regions_mouse <- fread("~/subprojects/chris_paper/phastcons/Mouse_enhancers_1kb_reformat.mm10.60way.phastCons.bed") %>% setnames(c("chr","start","end","name","class","conservation"))
all_regions_human <- fread("~/subprojects/chris_paper/phastcons/Human_ucsc_all_enhancers.hg38.60way.phastCons.bed") %>% setnames(c("chr","start","end","name","class","conservation"))
all_regions_mouse <- fread("~/subprojects/chris_paper/phastcons/Mouse_mm10_UCSC_all_enhancers.mm10.60way.phastCons.bed") %>% setnames(c("chr","start","end","name","class","conservation"))

outdir <- "~/subprojects/chris_paper/phastcons/"

if (species=="human"){
  regions <- regions_human
  regions_background <- all_regions_human
} else if (species=="mouse"){
  regions <- regions_mouse
  regions_background <- all_regions_mouse
}

# convert regions that are no in the phascons to 0 
regions[regions$conservation==".",]$conservation <- "0"
regions_background[regions_background$conservation==".",]$conservation <- "0"

# make conservation numeric
regions$conservation <- as.numeric(regions$conservation)
regions_background$conservation <- as.numeric(regions_background$conservation)

# split by primes and non promed and by tissue 
regions$group <- sapply(strsplit(regions$class,"_"),"[",2)
regions[regions$group=="prim",]$group <- "primed"
regions[regions$group=="nonprim",]$group <- "non-primed"
regions$tissue <- sapply(strsplit(regions$class,"_"),"[",1)

# add columns the background df 
regions_background$group <- "background"
regions_background$tissue <- "all"

# merge background and regions dfs
merged_df <- rbind(regions, regions_background)
merged_df$group <- factor(merged_df$group, levels=c("primed","non-primed","poised","background"))


# add the number of elements per group
df_summary <- as.data.frame(table(merged_df[,c("group")])) 
df_summary <- df_summary[df_summary$Freq!=0,]

stat <- ggpubr::compare_means(conservation ~ group,  data = merged_df, test="wilcox.test")
my_comparisons <- list( c("primed", "non-primed"), c("primed", "poised"), c("non-primed", "poised"), 
                          c("primed", "background"), c("poised", "background"), c("non-primed", "background") )

pdf(paste0(outdir,"/conservation_",species,".pdf"), width=3, height=4)
ggplot(merged_df, aes(x=group, y=conservation)) + 
  geom_violin() + geom_boxplot(width=0.1) + theme_bw() + ylab("Conservation Score") + xlab("") +  ylim(c(-0.05,2)) + ggtitle(species) + 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_text(data = df_summary, aes(label = Freq, y = -0.05), position = position_dodge(width = 0.9),size=2) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons, conservation ~ group,  data = merged_df, test="wilcox.test")
dev.off()

# add the number of elements per group
df_summary <- as.data.frame(table(merged_df[,c("group","tissue")])) 
df_summary <- df_summary[df_summary$Freq!=0,]

pdf(paste0(outdir,"/conservation_",species,"_split_by_tissue.pdf"), width=10, height=4)
ggplot(merged_df, aes(x=group, y=conservation, col=tissue)) + 
  geom_boxplot(position = position_dodge2(preserve = "single")) + theme_bw() + ylab("Conservation Score") + xlab("") + ylim(c(-0.05,1)) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  
  geom_text(data = df_summary, aes(label = Freq, y = -0.05), position = position_dodge(width = 0.9),size=1.75)
dev.off()

