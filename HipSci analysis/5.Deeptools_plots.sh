
module load python/3.9.7 
qlogin --mem=10Gb

grep "peaks_0_mut_1,peaks_1_mut_0" ../peaks_snvs_hipsci/interesting_motifs.txt  > mut_peaks_lost.txt
cut -f 13-23 mut_peaks_lost.txt | sed 's/peaks_0_mut_1/red/g' | sed 's/peaks_1_mut_1/black/g' | sed 's/peaks_1_mut_0/blue/g' | sed 's/peaks_0_mut_0/black/g'  > colours_lost.txt
cut -f 1 mut_peaks_lost.txt  | sed 's/_/\t/g' | awk '{print $1 "_" $2-1 "_" $2}' | sed 's/^/chr/g' > interesting_motifs_lost.bed
cut -f 3 mut_peaks_lost.txt  > which_enhancer_lost.bed
paste interesting_motifs_lost.bed which_enhancer_lost.bed colours_lost.txt > input_to_matrix_lost.txt


grep "peaks_0_mut_0,peaks_1_mut_1" ../peaks_snvs_hipsci/interesting_motifs.txt  > mut_peaks_gained.txt
cut -f 13-23 mut_peaks_gained.txt | sed 's/peaks_0_mut_1/black/g' | sed 's/peaks_1_mut_1/red/g' | sed 's/peaks_1_mut_0/black/g' | sed 's/peaks_0_mut_0/blue/g'  > colours_gained.txt
cut -f 1 mut_peaks_gained.txt  | sed 's/_/\t/g' | awk '{print $1 "_" $2-1 "_" $2}' | sed 's/^/chr/g' > interesting_motifs_gained.bed
cut -f 3 mut_peaks_gained.txt  > which_enhancer_gained.bed
paste interesting_motifs_gained.bed which_enhancer_gained.bed colours_gained.txt > input_to_matrix_gained.txt



# first annotate the bedgraphs with which enhancer they overlap
less input_to_matrix_lost.txt | cut -f 1  |  sed 's/_/\t/g' | sed 's/chr//g' | awk '{print $1 "\t" $2-2000 "\t" $3+2000}' > lost_sites.bed 

while read DON; do
    bedtools intersect -wa -wb -a lost_sites.bed -b ../bigwigs/${DON}_all_sum.bedgraph  > enhancer_overlap_${DON}_lost.txt
done < ../bigwigs/donors.txt 

# first annotate the bedgraphs with which enhancer they overlap
less input_to_matrix_gained.txt | cut -f 1  |  sed 's/_/\t/g' | sed 's/chr//g' | awk '{print $1 "\t" $2-2000 "\t" $3+2000}' > gained_sites.bed 

while read DON; do
    bedtools intersect -wa -wb -a gained_sites.bed -b ../bigwigs/${DON}_all_sum.bedgraph  > enhancer_overlap_${DON}_gained.txt
done < ../bigwigs/donors.txt 

mkdir -p plots

# make a matrix of accessibility
PHEN=$(echo "lost")
while read FILE; do 
    echo ${FILE} | cut -f 1 -d ' ' > file_name.tmp
    echo ${FILE} | cut -f 1 -d ' '  |  sed 's/_/\t/g' > tmp
    echo ${FILE} | cut -f 2 -d ' ' > enhancer_name.tmp
    ENH_ID=$(cut -f1,2,3 -d "_" enhancer_name.tmp)
    SNP_SITE=$(cut -f3 tmp)
    ENH_START=$(grep -w $ENH_ID ../enhancer_files/enhancer_files.txt | head -n1 | cut -f 2 | awk -v x=${SNP_SITE} '{print x-$1+1000}' | awk '{printf "%d00\n", $0 / 100}' )
    ENH_END=$(grep -w $ENH_ID ../enhancer_files/enhancer_files.txt | head -n1 | cut -f 3 | awk -v x=${SNP_SITE} '{print $1-x+1100}' |  awk '{printf "%d00\n", $0 / 100}' )
    echo ${FILE} | cut -f 3-12 -d ' ' |  sed 's/_/\t/g' > colours.tmp
    cat colours.tmp <(head -n1 ../peaks_snvs_hipsci/interesting_motifs.txt | cut -f 13-23 ) | awk '{ for (i=1; i<=NF; i++) RtoC[i]= (i in RtoC?RtoC[i] OFS :"") $i; }  END{ for (i=1; i<=NF; i++) print RtoC[i] }' | grep red | cut -f 2 -d " "   > mutation_donors.txt 
    cat colours.tmp <(head -n1 ../peaks_snvs_hipsci/interesting_motifs.txt | cut -f 13-23 ) | awk '{ for (i=1; i<=NF; i++) RtoC[i]= (i in RtoC?RtoC[i] OFS :"") $i; }  END{ for (i=1; i<=NF; i++) print RtoC[i] }' | grep blue | cut -f 2 -d " "  > no_mutation_donors.txt 
    while read RED_DON; do 
        echo ${RED_DON}
        CHR=$(cut -f 1 tmp | sed 's/chr//g')
        START=$( awk '{print $2-2000}' tmp)
        END=$( awk '{print $3+2000}' tmp)
        awk -v x=${RED_DON} -v c=${CHR} -v p=${START} -v q=${END} '$1==c && $2==p && $3==q {print x "\t" $5 "\t" $6 "\t" $7}' enhancer_overlap_${RED_DON}_${PHEN}.txt >> mutation_donors_counts.txt
    done < mutation_donors.txt
    while read BLUE_DON; do 
        echo ${BLUE_DON}
        CHR=$(cut -f 1 tmp | sed 's/chr//g')
        START=$( awk '{print $2-2000}' tmp)
        END=$( awk '{print $3+2000}' tmp)
        awk -v x=${BLUE_DON} -v c=${CHR} -v p=${START} -v q=${END} '$1==c && $2==p && $3==q {print x "\t" $5 "\t" $6 "\t" $7}' enhancer_overlap_${BLUE_DON}_${PHEN}.txt >> no_mutation_donors_counts.txt
    done < no_mutation_donors.txt
    sed 's/Bima1/1/g' mutation_donors_counts.txt | sed 's/Kolf2/2/g' | sed 's/Kolf3/3/g'  | sed 's/Kucg2/4/g'  | sed 's/Letw5/5/g' | sed 's/Podx1/6/g' | sed 's/Qolg1/7/g'  | sed 's/Sojd3/8/g' | sed 's/Wibj2/9/g' | sed 's/Yoch6/10/g' | cat - chr_for_overlap.txt | sort -k1,1 -k2,2n  > mutation_donors_artifical_coord.bed
    sed 's/Bima1/1/g' no_mutation_donors_counts.txt | sed 's/Kolf2/2/g' | sed 's/Kolf3/3/g'  | sed 's/Kucg2/4/g'  | sed 's/Letw5/5/g' | sed 's/Podx1/6/g' | sed 's/Qolg1/7/g'  | sed 's/Sojd3/8/g' | sed 's/Wibj2/9/g' | sed 's/Yoch6/10/g' | cat - chr_for_overlap.txt | sort -k1,1 -k2,2n > no_mutation_donors_artifical_coord.bed
    bedGraphToBigWig mutation_donors_artifical_coord.bed artifical_coord_chrom.sizes mutation_donors_artifical_coord.bw
    bedGraphToBigWig no_mutation_donors_artifical_coord.bed artifical_coord_chrom.sizes no_mutation_donors_artifical_coord.bw
    MAX_COUNT=$( cat mutation_donors_artifical_coord.bed no_mutation_donors_artifical_coord.bed | sort -nrk4,4 | head -n1 | awk '{print $4+10}' )
    awk -v x=$(cut -f 2 tmp) -v y=$(cut -f 3 tmp) '{print $1 "\t" x "\t" y }' <( cut -f 1 chr_for_overlap.txt) > regions.tmp 
    FILE_NAME=$(less file_name.tmp)
    NO_MUT=$( wc -l no_mutation_donors.txt | cut -f 1 -d " ")
    MUT=$( wc -l mutation_donors.txt| cut -f 1 -d " ")
    computeMatrix reference-point -R regions.tmp -o ${FILE_NAME}.matrix -b ${ENH_START} -a ${ENH_END} \
                --binSize 50 \
                --skipZeros \
                -S mutation_donors_artifical_coord.bw no_mutation_donors_artifical_coord.bw \
                --samplesLabel "mut ($MUT)" "ref ($NO_MUT)"
    ENH_NAME=$(cat enhancer_name.tmp)
    plotProfile --matrixFile ${FILE_NAME}.matrix    \
            -out plots_sd/${FILE_NAME}_${PHEN}_merged.pdf \
            --colors red blue --perGroup \
            --legendLocation best --plotHeight 10 --plotHeight 10 --plotTitle ${ENH_NAME}_${PHEN} \
            --startLabel enhStart --endLabel enhEnd --averageType mean\
            --plotType std --yMin 0 0 --yMax ${MAX_COUNT} ${MAX_COUNT}
    rm -rf mutation_donors_counts.txt no_mutation_donors_counts.txt 
    rm -rf *tmp* mutation_donors* no_mutation_donors*
done < input_to_matrix_lost.txt

























