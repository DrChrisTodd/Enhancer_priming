
export PATH=/bi/group/reik/Jannat/nmt/analysis/methacc/meth/motif_analysis/meme-5.5.0/bin:/bi/group/reik/Jannat/nmt/analysis/methacc/meth/motif_analysis/meme-5.5.0/libexec/meme-5.5.0:$PATH

module load bedtools 
bedtools getfasta -fi genome/Mus_musculus_GRCm38.fa -bed <( sed 's/chr//g' ../methylation_states_at_enhancers/refilter/E7.5_Ect_Enh_edit.refiltered2.bed) > E7.5_Ect_Enh_edit.refiltered2.fa 
bedtools getfasta -fi genome/Mus_musculus_GRCm38.fa -bed <( sed 's/chr//g' ../../../H3K27ac/overlapping_enhancers/adult_brain_enhancers_overlap.bed) > adult_brain_enhancers_overlap.bed 

./../../../../nmt/analysis/methacc/meth/motif_analysis/meme-5.5.0/bin/sea --p ../homer_calls/motifs_new_sequences.fa --m ../../../../nmt/analysis/methacc/meth/motif_analysis/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme 
./../../../../nmt/analysis/methacc/meth/motif_analysis/meme-5.5.0/bin/sea --p ../homer_calls/no_motifs_new_sequences.fa --m ../../../../nmt/analysis/methacc/meth/motif_analysis/motif_databases/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme  -n motifs_old_sequences.fa
./../../../../nmt/analysis/methacc/meth/motif_analysis/meme-5.5.0/bin/sea --p ../homer_calls/motifs_new_sequences.fa --m ../../../../nmt/analysis/methacc/meth/motif_analysis/motif_databases/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme -n no_motifs_old_sequences.fa

awk '{print $11 "\t" $4}' ../motifs_new_sequences.txt | sed 's/^/>/g' | sed 's/\t/\n/g' > motifs_old_sequences.fa 

findMotifs.pl -find knownResults/known1.motif motifs_new_sequences.fa fasta > initial_motifs_altered_known.txt

while read FILE; do 
    findMotifs.pl motifs_new_sequences.fa human initial_motifs_altered_rerun/ -find initial_motifs_altered/knownResults/${FILE}.motif > initial_motifs_altered/initial_motifs_altered_${FILE}.txt
done < known_motifs.txt
