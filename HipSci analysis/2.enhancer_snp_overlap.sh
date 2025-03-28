while read DONOR; do 
    bedtools intersect -wa -wb -a <( awk -v donor=${DONOR} '{print $1 "\t" $2 "\t" $2 "\t" donor}' ../HipSci/Variants/${DONOR}.vcf | grep -v ^# ) -b   enhancer_atac_overlap/enhancers_with_at_least1overlap_coord.txt >  ${DONOR}_enhancer_overlap.txt
done <donors.txt