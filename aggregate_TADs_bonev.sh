# Aggregate plots of TADs from Bonev et al. 2017 (https://pubmed.ncbi.nlm.nih.gov/29053968/)
# 1) Lift TADs from mm10 to mm39 using liftOver and GenomicRanges (in R)
# 2) Split autosomes and chrX
# 3) Save in files Bonev_ESC_TADs_mm39_AUTOSOMES.bed and Bonev_ESC_TADs_mm39_chrX.bed

cd working_dir
mkdir -p aggregate_plots

module load Python/3.9.6-GCCcore-11.2.0

TADS_AUTOSOMES=path_to/Bonev_ESC_TADs_mm39_AUTOSOMES.bed
TADS_X=path_to/Bonev_ESC_TADs_mm39_chrX.bed


for sample in female male ; do
  for cond in GFPpos GFPneg ; do

    hic_matrix=path_to_hic_matrices/merged_$sample\_$cond\_10kb.hic 

    fanc aggregate $hic_matrix $TADS_AUTOSOMES \
    aggregate_plots/merged_$sample\_$cond\_10kb_AUTOSOMES.agg \
    -p aggregate_plots/merged_$sample\_$cond\_10kb_AUTOSOMES.pdf \
    -m aggregate_plots/merged_$sample\_$cond\_10kb_AUTOSOMES.agg.txt \
    --tad-strength aggregate_plots/merged_$sample\_$cond\_10kb_AUTOSOMES.TADstrength.txt \
    -e -r 1.0 --rescale --vmin 0.03 --vmax 0.125

    fanc aggregate $hic_matrix $TADS_X \
    aggregate_plots/merged_merged_$sample\_$cond\_10kb_chrX.agg \
    -p aggregate_plots/merged_$sample\_$cond\_10kb_chrX.pdf \
    -m aggregate_plots/merged_$sample\_$cond\_10kb_chrX.agg.txt \
    --tad-strength aggregate_plots/merged_$sample\_$cond\_10kb_chrX.TADstrength.txt \
    -e -r 1.0 --rescale --vmin 0.03 --vmax 0.125

  done
done




