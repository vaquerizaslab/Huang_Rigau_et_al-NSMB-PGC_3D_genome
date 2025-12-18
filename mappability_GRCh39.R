# Low mappability regions

library(GenomicRanges)

# File k36.merged.GRCm39_umap.bed shows uniquely mapping regions. 
# These were downloaded from https://bismap.hoffmanlab.org/ (Karimzadeh et al., 2018) 


# UMAP bed file (uniquely mapping loci)
k36 <- read.table("data/k36.merged.GRCm39_umap.bed", skip = 1)
k36 <- makeGRangesFromDataFrame(k36, seqnames.field = "V1", start.field = "V2", end.field = "V3")

chrom_sizes <- read.table("data/chromSizes_GRCm39.tsv", row.names = 1)
sizes <- chrom_sizes$V3
names(sizes) <- rownames(chrom_sizes)


# Calculate how much of a 100kb window is uniquely mapping
for(binsize in c( 100000)) {
  
  if(binsize == 100000) simpl_binsize <- "100kb"
  # if(binsize == 500000) simpl_binsize <- "500kb"
  # if(binsize == 1000000) simpl_binsize <- "1Mb"
  
  genome <- tileGenome(sizes, tilewidth = binsize, cut.last.tile.in.chrom = T)
  genome$percent_uniq_map <- as.numeric(NA)
  # I also add this column to k36, otherwise I cannot append the two objects
  k36$percent_uniq_map <- NA
  
  k36_low_map <- GRanges() # List of regions with less than 90% uniquely mapping coverage
  
  for(chr in paste0("chr", c(1:19, "X","Y"))) {
    chr_bins <- genome[seqnames(genome) == chr]
    uniq_map <- c()

    # Here I append 1 bin to uniq_map calculate the coverage of the new object.
    # For separate objects coverage is maximum 1, so it will be 2 only in the regions from the bin that have uniquely mapping regions.
    # Thus, I sum the lengths that have coverage = 2.

    for(i in 1:length(chr_bins)) {
      uniq_map <- c(uniq_map,
                    suppressWarnings(sum(coverage(append(chr_bins[i],
                                                         k36))[[chr]]@lengths[coverage(append(chr_bins[i],
                                                                                              k36))[[chr]]@values == 2])))

      message(paste(chr, " - bin", i ,"of", length(chr_bins)))
    }
    genome[seqnames(genome) == chr]$percent_uniq_map <- uniq_map/binsize*100
    message(chr)

    per_chr <- genome[seqnames(genome) == chr]

    # Get regions with less than 90% unique sequence.
    # Seems very strict, but was useful to filter out bad loop calls and keep good ones, after testing other thresholds
    low_map_per_chr <- per_chr[per_chr$percent_uniq_map < 90]
    k36_low_map     <- append(k36_low_map, low_map_per_chr)
    
  }
  

  # Save RData of uniquely mapping regions 
  # save(genome, file = paste0("output/k36_genome_percent_uniquely_mapping_per_", simpl_binsize, "_bins.RData"))
  write.table(as.data.frame(k36_low_map)[, c(1,2,3)],
              file = paste0("processed_data/k36_", simpl_binsize, "_regions_less_90percent_unique.bed"),
              sep = "\t", quote = F, col.names = F, row.names = F)
}


