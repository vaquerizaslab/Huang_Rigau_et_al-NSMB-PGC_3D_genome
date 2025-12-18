# Loop call filtering

# Loop calls were obtained as follows
# 1) Merged replicate matrices were transformed from .hic to .mcool using fanc to-cooler
# 2) Mustache was run on these matrices with the following parameters
#     --sparsityThreshold 0.9 \
#     --resolution 10kb \
#     --pThreshold 0.01 \
#     --sigmaZero 1.6 \
#     -p 1 



library(GenomicRanges)

blacklist <- makeGRangesFromDataFrame(read.table("data/GRCm39-blacklist_lifted_from_mm10_merged.bed"),
                                      seqnames.field = "V1", start.field = "V2", end.field = "V3")

# Regions with low mappability based on script mappability_GRCh39.R
k36_low_map <- makeGRangesFromDataFrame(read.table("processed_data/k36_100kb_regions_less_90percent_unique.bed"),
                                        seqnames.field = "V1", start.field = "V2", end.field = "V3")
# Unify both low mappability lists.
low_map <- reduce(append(blacklist, k36_low_map))


# I need to select loops in autosomes
# Remove the ones in blacklisted regions or repetitive regions
# Save score so I can visualize it in FANC
# For selecting windows to plot, resize to 1Mb unless they were bigger! 


for(samplename in c("female_GFPpos", "female_GFPneg", "male_GFPpos", "male_GFPneg")) {

  loop_calls <- read.table(paste0("data/raw_loops/cooler_out_", samplename, ".tsv"), header = T)

  
  # Checking that none of the anchors overlaps with blacklisted regions
  anchor1 <- makeGRangesFromDataFrame(loop_calls, seqnames.field = "BIN1_CHR", start.field = "BIN1_START", end.field = "BIN1_END")
  anchor2 <- makeGRangesFromDataFrame(loop_calls, seqnames.field = "BIN2_CHROMOSOME", start.field = "BIN2_START", end.field = "BIN2_END")
  
  # Low mappability filtering
  loops_filt <- loop_calls[!overlapsAny(anchor1, low_map) & !overlapsAny(anchor2,low_map), ]
  
  # Autosomes only
  loops_filt2 <- loops_filt[loops_filt$BIN1_CHR %in% paste0("chr", 1:19), ]
  


  # CENTER TO CENTER BED FILE
  ctoc <- data.frame(chr = loops_filt2$BIN1_CHR,
                     start = apply(loops_filt2[, c("BIN1_START", "BIN1_END")], 1, function(x) round(mean(x))),
                     end   = apply(loops_filt2[, c("BIN2_START", "BIN2_END")], 1, function(x) round(mean(x))),
                     name = paste0("loop_", 1:nrow(loops_filt2)), 
                     score = loops_filt2$DETECTION_SCALE,
                     strand = rep("+", nrow(loops_filt2)))
  
  
  write.table(ctoc,
              col.names = F, row.names = F, quote = F, sep = "\t",
              file = paste0("processed_data/filtered_loops_in_", samplename, "_center_to_center_AUTOSOMES_noBLACKLIST_noLOWMAP.bed"))
  write.table(loops_filt2,
              col.names = F, row.names = F, quote = F, sep = "\t",
              file = paste0("processed_data/filtered_loops_in_", samplename, "_AUTOSOMES_noBLACKLIST_noLOWMAP.bedpe"))
  
}


