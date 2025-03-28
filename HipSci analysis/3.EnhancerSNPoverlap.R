# Rscript for decipfer enhancers - for chris


library(data.table)
library(purrr)
library(ggplot2)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)

mypaste <- function(x,y) paste(x, y, sep="_")


io <- list()
#output directory
io$outdir <- "/bi/group/reik/Jannat/chris_snps_atac/peaks_snvs_hipsci/"
io$enhancer_snv_list <- "/bi/group/reik/Jannat/chris_snps_atac/peaks_snvs_hipsci/enhancer_snv_overlap/enhancers_overlapping_snps.txt"
io$atac_enhancer_overlap_dir <-  "/bi/group/reik/Jannat/chris_snps_atac/peaks_snvs_hipsci/enhancer_atac_overlap/"
io$variant_enhancer_overlap_dir <-  "/bi/group/reik/Jannat/chris_snps_atac/peaks_snvs_hipsci/enhancer_snv_overlap/"
#io$motifs_collapsed <-  "/bi/group/reik/Jannat/chris_snps_atac/peaks_snvs_hipsci/snp_motif_overlap/snp_motif_overlap_collapsed.txt"
#io$motifs <-  "/bi/group/reik/Jannat/chris_snps_atac/peaks_snvs_hipsci/snp_motif_overlap/snp_motif_overlap.txt"

opts <- list()
# list of enhancers whihc overlap atac peaks files
opts$files_actac_overlap <- list.files(io$atac_enhancer_overlap_dir,pattern="*_with_overlapped_region.bed")
# list of enhancers whihc overlap atac snps files
opts$files_snv_overlap <- list.files(io$variant_enhancer_overlap_dir,pattern="HPSI")

# read in the atac overlaps with enhancers 
file_atac <- lapply(opts$files_actac_overlap, function(file) fread(sprintf("%s/%s",io$atac_enhancer_overlap_dir,file), stringsAsFactors=T)) %>%  
    rbindlist %>%  setnames(c("chr_peak","start_peak","end_peak","peak_id","info","strand","val1","val2","val3","val4","chr_enh","start_enh","end_enh","enh_id"))
file_atac$donor_sample = paste(as.character(lapply(strsplit(as.character(file_atac$peak_id), split="_"), "[", c(3))), as.character(lapply(strsplit(as.character(file_atac$peak_id), split="_"), "[", c(4))), sep="_")
file_atac$donor = as.character(lapply(strsplit(as.character(file_atac$peak_id), split="_"), "[", c(3)))

# read in the snv overlaps with enhancers 
file_snv <- lapply(opts$files_snv_overlap, function(file) fread(sprintf("%s/%s",io$variant_enhancer_overlap_dir,file), stringsAsFactors=T)) %>%  
    rbindlist %>%  setnames(c("chr_snp","start_snp","end_snp","var_id","chr_enh","start_enh","end_enh","enh_id"))
file_snv$donor <- str_to_title(paste(substr(as.character(file_snv$var_id),11,14),substr(as.character(file_snv$var_id),16,16),sep=""))
file_snv$snp_id <- paste(file_snv$chr_snp, file_snv$start_snp,sep="_")

# read in a list of snps and enhancer overlaps
enhancer_snv_list<-c(fread(io$enhancer_snv_list,stringsAsFactors=F,header=F)$V1)
#enhancer_snv_list<-enhancer_snv_list[1:20]
#get a list of the snps 
snv<- unlist(strsplit(enhancer_snv_list,split="-"))[2*(1:length(enhancer_snv_list))-1]
#get a list of the enhancer
enhancers<- unlist(strsplit(enhancer_snv_list,split="-"))[2*(1:length(enhancer_snv_list))  ]

#get a list of the donor ids
donors <- sort(unique(file_atac$donor))

#make a df to populate with whether there is a peak per donor
df <- data.table(enhancer_snv_list,snv,enhancers) %>% setnames(c("id","snv","enhancer")) %>% .[,as.vector(donors):=as.character("0")]

# round 1 look through atac peaks to decide if there is a peak there 
for (enh in enhancers) {
    #enh_tmp = paste(strsplit(enh,split="_")[[1]][1:length( strsplit(enh,split="_")[[1]])-1],collapse="_") # remove the annotation from the enhancer 
    # get which atac peaks overlap the enhancer
    atac_peaks_for_enh <- unique(file_atac[file_atac$enh==enh,c("enh_id","donor_sample","donor")])
    # get a list of which donors get a peak 
    donors_with_peak <- atac_peaks_for_enh$donor
    # add all the donors to the list so that each donor appears once in the counts table
    donors_with_peak <- c(donors_with_peak,donors)
    # get a count table
    counts<-table(donors_with_peak)
    #for (d in names(counts)){
    #    df[enhancer==enh,c(d):=paste("peaks_",as.character(counts[[d]]-1),sep="")]
    #}
    for (d in names(counts)){
        # if the donor appears more than once then there is a peak it is present three then there are peaks in the at least 2 out of 3 of the repeats of samples
        if (counts[[d]]>=3){ 
            df[enhancer==enh,c(d):="peaks_1"]
        # if the donor appears once then there is not a peak 
        } else {
             df[enhancer==enh,c(d):="peaks_0"]
        }
    }
}

fwrite(df, file=paste0(io$outdir,"/enh.tmp"), sep="\t", row.names=F, col.names=T)
df <- fread(paste0(io$outdir,"/enh.tmp"), stringsAsFactors=FALSE)


#make a df to populate with whether there is a snp in the donor
df2 <- data.table(enhancer_snv_list,snv,enhancers) %>% setnames(c("id","snv","enhancer")) %>% .[,as.vector(donors):=as.character("0")]

for (mut in snv) {
    # get which snps are found in the enhancer
    snv_peaks_for_enh <- file_snv[file_snv$snp_id==mut]
    # get a list of which donors have an sv  
    donors_with_peak <- snv_peaks_for_enh$donor
    # add all the donors to the list so that each donor appears once in the counts table
    donors_with_peak <- c(donors_with_peak,donors)
    # get a count table
    counts<-table(donors_with_peak)
    #for (d in names(counts)){
    #    df2[snv==mut,c(d):=paste("mut_",as.character(counts[[d]]-1),sep="")]
    #}
    for (d in names(counts)){
        if (counts[[d]]>1){
        # if the donor appears more than once then there is a snv
            df2[snv==mut,c(d):="mut_1"]
        # if the donor appears once then there is not a snv
        } else {
            df2[snv==mut,c(d):="mut_0"]
        }
    }
}

fwrite(df2, file=paste0(io$outdir,"/mut.tmp"), sep="\t", row.names=F, col.names=T)
df2 <- fread(paste0(io$outdir,"/mut.tmp"), stringsAsFactors=FALSE)


# merge the atac df and the snv df so each donor column says whether there is a peak and a snp
final_df<- cbind(df[,1:3],mapply(mypaste, df[,4:13],df2[,4:13]))

# wirte df out 
fwrite(final_df, file=paste0(io$outdir,"/mut_peaks_df_binary.txt"), sep="\t", row.names=F, col.names=T)

final_df <- fread(paste0(io$outdir,"/mut_peaks_df_binary.txt"), stringsAsFactors=FALSE)




## remove the donor information and just get the counts of the different peak snv combinations
# make a df that you can populate with the cound of each peak snv combination 
unique_df <- data.table(enhancer_snv_list,snv,enhancers) %>% setnames(c("id","snv","enhancer")) %>% .[,"unique":=as.character("0")] %>% .[,c("peaks_0_mut_0","peaks_0_mut_1","peaks_1_mut_0","peaks_1_mut_1"):=as.numeric(0)]
for (row in 1:nrow(unique_df)){
    # sort the order of the peak mut combinations 
    peaks_muts=sort(as.character(final_df[row,4:13]))
    # get all possible peak snv combinations
    unique_df[row,4]=paste(unique(peaks_muts),collapse=",")
    # get counts of peak snv combinations
    tabled<-table(peaks_muts)
    # fill in the df
    for (item in unique(peaks_muts)){
        unique_df[row,c(item):=tabled[[item]]] 
    }
    #count_table=c()
    #for (item in 1:dim(tabled)){
    #    count_table=c(count_table,paste0(names(tabled)[item],":",unname(tabled)[item]))
    #}
    #unique_df[row,5]=paste(unique(count_table),collapse=",")
}

fwrite(unique_df, file=paste0(io$outdir,"/mut_peaks_df_unique_binary.txt"), sep="\t", row.names=F, col.names=T)

# addd enhancer group information
unique_df$group <- unlist(strsplit(unique_df$enhancer,split="_"))[4*(1:length(unique_df$enhancer))  ]
# get the number of different peak snv combinations
unique_df$length<-str_count(unique_df$unique,",")+1
# remove H1 and MSC enhancers
unique_df <- unique_df[unique_df$group!="H1" & unique_df$group!="MSC"]

# allow some error in the mutation calling 
# generate PASS or fail based on consistency of the peak snv consistency
unique_df$QC<-"NA"
for (row in 1:nrow(unique_df)){
    if (unique_df$length[row]==4 | unique_df$length[row]==1){
        unique_df$QC[row]<-"FAIL"
    } else if (unique_df$length[row]==2){
        unique_df$QC[row]<-"PASS"
    } else {
        if (unique_df$peaks_0_mut_0[row]==1 & unique_df$peaks_0_mut_1[row]!=1 & unique_df$peaks_1_mut_0[row]!=1 & unique_df$peaks_1_mut_1[row]!=1){
            unique_df$QC[row]<-"PASS"
            unique_df$unique[row]<- gsub("peaks_0_mut_0,", "",  unique_df$unique[row])
            unique_df$unique[row]<- gsub(",peaks_0_mut_0", "",  unique_df$unique[row])
        } else if (unique_df$peaks_0_mut_0[row]!=1 & unique_df$peaks_0_mut_1[row]==1 & unique_df$peaks_1_mut_0[row]!=1 & unique_df$peaks_1_mut_1[row]!=1){
            unique_df$QC[row]<-"PASS"
            unique_df$unique[row]<- gsub("peaks_0_mut_1,", "",  unique_df$unique[row])
            unique_df$unique[row]<- gsub(",peaks_0_mut_1", "",  unique_df$unique[row])
        } else if (unique_df$peaks_0_mut_0[row]!=1 & unique_df$peaks_0_mut_1[row]!=1 & unique_df$peaks_1_mut_0[row]==1 & unique_df$peaks_1_mut_1[row]!=1){
            unique_df$QC[row]<-"PASS"
            unique_df$unique[row]<- gsub("peaks_1_mut_0,", "",  unique_df$unique[row])
            unique_df$unique[row]<- gsub(",peaks_1_mut_0", "",  unique_df$unique[row])
        } else if (unique_df$peaks_0_mut_0[row]!=1 & unique_df$peaks_0_mut_1[row]!=1 & unique_df$peaks_1_mut_0[row]!=1 & unique_df$peaks_1_mut_1[row]==1){
            unique_df$QC[row]<-"PASS"
            unique_df$unique[row]<- gsub("peaks_1_mut_1,", "",  unique_df$unique[row])
            unique_df$unique[row]<- gsub(",peaks_1_mut_1", "",  unique_df$unique[row])
        } else {
            unique_df$QC[row]<-"FAIL"
        }
    }
}

unique_df$confidence<-0
# anything above 0.4 is confident - ie 2/10 supports for the read 
for (row in 1:nrow(unique_df)){
    if (unique_df$length[row]==4 | unique_df$length[row]==1){
        unique_df$confidence[row]<-0
    } else if (unique_df$length[row]==2){
        annos<-unlist(strsplit(unique_df$unique[row],split=","))
        count1<-as.numeric(unique_df %>% select(annos[1]) %>% .[row,])
        count2<-as.numeric(unique_df %>% select(annos[2]) %>% .[row,])
        if(count1<count2) {
            unique_df$confidence[row]<-(count1/(count1+count2))*2
        } else {
            unique_df$confidence[row]<-(count2/(count1+count2))*2
        }
    } else if (unique_df$length[row]==3 & unique_df$QC[row]=="PASS") {
        annos<-unlist(strsplit(unique_df$unique[row],split=","))
        count1<-as.numeric(unique_df %>% select(annos[1]) %>% .[row,])
        count2<-as.numeric(unique_df %>% select(annos[2]) %>% .[row,])
        if(count1<count2) {
            unique_df$confidence[row]<-(count1/(count1+count2))*2
        } else {
            unique_df$confidence[row]<-(count2/(count1+count2))*2
        }
     } else if (unique_df$length[row]==3 & unique_df$QC[row]=="FAIL"){
            unique_df$confidence[row]<-0
    } else {
        print("you missed something")
    }
}

fwrite(unique_df, file=paste0(io$outdir,"/mut_peaks_df_unique_binary_full_info.txt"), sep="\t", row.names=F, col.names=T)

