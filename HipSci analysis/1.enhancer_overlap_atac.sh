awk '$5!="Brain"' ../../enhancer_files/enhancer_files.txt | awk '$5!="NPC"' | awk '$5!="ME"'  > enhancers_unique.bed
while read FILE; do bedtools intersect -wa -wb -a enhancers_unique.bed -b ../../HipSci/ATAC_peaks/${FILE}_GRCh38_bowtie2_peaks_no_mt.narrowPeak > ${FILE}_enhancer_peak_overlap.bed; done < peak_files.txt 
cat *overlap.bed | cut -f 4 | sort | uniq > enhancers_with_at_least1overlap.txt
cat *overlap.bed | cut -f 1,2,3,4,5 | sort | uniq > enhancers_with_at_least1overlap_coord.txt