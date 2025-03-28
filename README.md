# Enhancer_priming

Repository for scripts used in the bioinformatic analyses of the following publication:

Functional evaluation of transposable elements as enhancers in mouse embryonic and trophoblast stem cells

Christopher D Todd, Jannat Ijaz, Fereshteh Torabi, Oleksandr Dovgusha, Stephen Bevan, Olivia Cracknell, Tim Lohoff, Stephen Clark, Ricard Argelaguet, Juliette Pierce, Ioannis Kafetzopoulos, Alice Santambrogio, Jennifer Nichols, Ferdinand von Meyenn1, Ufuk Günesdogan, Stefan Schoenfelder, Wolf Reik

doi:  https://doi.org/10.1101/2024.09.09.611867 


## Abstract

Embryonic development requires the accurate spatiotemporal execution of cell lineage-specific gene expression programs, which are controlled by transcriptional enhancers. Developmental enhancers adopt a primed chromatin state prior to their activation; however how this primed enhancer state is established, maintained, and how it affects the regulation of developmental gene networks remains poorly understood. Here, we use comparative multi-omic analyses of human and mouse early embryonic development to identify subsets of post-gastrulation lineage-specific enhancers which are epigenetically primed ahead of their activation, marked by the histone modification H3K4me1 within the epiblast. We show that epigenetic priming occurs at lineage-specific enhancers for all three germ layers, and that epigenetic priming of enhancers confers lineage-specific regulation of key developmental gene networks. Surprisingly in some cases, lineage-specific enhancers are epigenetically marked already in the zygote, weeks before their activation during lineage specification. Moreover, we outline a generalisable strategy to use naturally occurring human genetic variation to delineate important sequence determinants of primed enhancer function. Our findings identify an evolutionarily conserved program of enhancer priming and begin to dissect the temporal dynamics and mechanisms of its establishment and maintenance during early mammalian development.


### Notes
Coordiantes for epiblast primed (ePrimed), epiblast poised (ePoised), and epiblast non-primed (eNonprimed) enhancers are in "./Primed_enhancer_annotations/"
Human coordinates are for hg38, and mouse are for mm10. 

The following scripts were used in the generation of the main figures:
- 1_Defining_primed_enhancers.Rmd
- 2_Homer_ChIP-seq_and_ATAC_analysis.Rmd
- 3_DNA_methylation_Acc_statistical_tests.R
- 4_PCHiC_gene_calling.Rmd
- 5_Plotting_gene_expression.R
- 6_Gene_ontology_enrichment_enrichr.R

Analysis of SNPs within HipSci iPSC lines are contained within "./HipSci analysis/":
1. enhancer_overlap_atac.sh - generated which atac peaks are in which donors
2. enhancer_snp_overlap.sh - generated which snps are in which donors
3. EnhancerSNPoverlap.R - generates sites which overlap atac peaks and snps which can be used to identify interesting sites
4. DecipeherEnhancerOverlaps.R - fishers exact of enrichment of the changes in specific motifs
5. Deeptools_plots.sh - generate example sites
6. findmotifs.sh - generate motifs at sites of interest

The following scripts were used in the generation of supplimentary figures/revision analyses:
- Getting_poised_enhancers.R
- Alluvial groupings.R
- Poised enhancer sankeys.R
- Enhancer decomissioning ChIPseq counts.R
- TF network analysis.R
- Analysis of enhancer sequence conservation phastcon.R
