# An insulation bin size of 300kb was selected for boundary calling after visual inspection of other bin sizes (from 70kb to 500 kb)

library(GenomicRanges)
library(gtools)
library(ggplot2)
library(devtools)


female_soma <- makeGRangesFromDataFrame(read.table("data/raw_boundaries/merged_female_GFPneg_10kb.insulation_boundaries_300kb_window.bed", col.names = c("chr", "start", "end", "x", "score", "strand")), keep.extra.columns = T)
female_PGC  <- makeGRangesFromDataFrame(read.table("data/raw_boundaries/merged_female_GFPpos_10kb.insulation_boundaries_300kb_window.bed", col.names = c("chr", "start", "end", "x", "score", "strand")), keep.extra.columns = T)
male_soma   <- makeGRangesFromDataFrame(read.table("data/raw_boundaries/merged_male_GFPneg_10kb.insulation_boundaries_300kb_window.bed",   col.names = c("chr", "start", "end", "x", "score", "strand")), keep.extra.columns = T)
male_PGC    <- makeGRangesFromDataFrame(read.table("data/raw_boundaries/merged_male_GFPpos_10kb.insulation_boundaries_300kb_window.bed",   col.names = c("chr", "start", "end", "x", "score", "strand")), keep.extra.columns = T)

# Low mappability 
# wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz
# Lifted to mm39 using http://genome-euro.ucsc.edu/cgi-bin/hgLiftOver
low_map <- makeGRangesFromDataFrame(read.table("data/GRCm39-blacklist_lifted_from_mm10_merged.bed", col.names = c("chr", "start", "end")), keep.extra.columns = T)
# I resize low-mappability regions (add 100kb, 50kb upstream and 50kb downstream)
start(low_map) <- start(low_map) - 50000
end(low_map) <- end(low_map) + 50000


### OBTAINING A FINAL LIST OF BOUNDARIES 

table_types_of_boundaries <- c()
for(sex in c("male", "female")) {
  boundaries_soma <- get(paste0(sex, "_soma"))
  boundaries_PGC <- get(paste0(sex, "_PGC"))
  
  # Make boundaries bigger so that I catch overlap with nearby boundaries
  start(boundaries_soma) <- start(boundaries_soma) - 15000
  end(boundaries_soma) <- end(boundaries_soma) + 15000
  
  start(boundaries_PGC) <- start(boundaries_PGC) - 15000
  end(boundaries_PGC) <- end(boundaries_PGC) + 15000
  
  # Append and merge
  all_boundaries <- append(boundaries_soma, boundaries_PGC)
  all_boundaries_merged <- reduce(all_boundaries)
  
  # Scores (default score is 0)
  all_boundaries_merged$score_in_PGC  <- 0
  all_boundaries_merged$score_in_soma <- 0
  
  # Calculate score per sample
  ## 1. calculate overlap
  fo_soma <- findOverlaps(all_boundaries_merged, boundaries_soma)
  ## 2. Get maximum score in case more than one TAD boundary overlaps
  max_score <- by(boundaries_soma$score[ subjectHits(fo_soma)], queryHits(fo_soma), max)
  ## 3. assign score to all merged object
  all_boundaries_merged$score_in_soma[as.integer(names(max_score))] <- unlist(max_score)
  
  # Repeat for PGCs 
  fo_PGC <- findOverlaps(all_boundaries_merged, boundaries_PGC)
  max_score <- by(boundaries_PGC$score[ subjectHits(fo_PGC)], queryHits(fo_PGC), max)
  all_boundaries_merged$score_in_PGC[as.integer(names(max_score))] <- unlist(max_score)
  
  # Remove lines with NA
  all_boundaries_merged <- all_boundaries_merged[!(is.na(all_boundaries_merged$score_in_PGC) |
                                                     is.na(all_boundaries_merged$score_in_soma))]
  
  # Select only boundaries with score >= 0.3 in at least one sample
  all_boundaries_merged <- all_boundaries_merged[all_boundaries_merged$score_in_PGC >= 0.3 | all_boundaries_merged$score_in_soma >= 0.3]
  
  # Autosomes only
  all_boundaries_merged <- all_boundaries_merged[seqnames(all_boundaries_merged) %in% paste0("chr", 1:19)]
  
  # Remove boundaries overlapping low-mappability regions 
  all_boundaries_merged <- all_boundaries_merged[!overlapsAny(all_boundaries_merged, low_map)]
  
  # Reduce size
  start(all_boundaries_merged) <- start(all_boundaries_merged) + 15000
  end(all_boundaries_merged)   <- end(all_boundaries_merged) - 15000
  
  # Classification 1
  shared <- all_boundaries_merged[all_boundaries_merged$score_in_PGC > 0 & all_boundaries_merged$score_in_soma > 0]
  pgc    <- all_boundaries_merged[all_boundaries_merged$score_in_PGC > 0 & all_boundaries_merged$score_in_soma == 0]
  soma   <- all_boundaries_merged[all_boundaries_merged$score_in_PGC == 0 & all_boundaries_merged$score_in_soma > 0]
  
  
  ## ADD in shared group if they are stronger in soma or in PGC  
  # - PGC-specific
  # - PGC-stronger
  # - Soma-stronger
  # - Soma-specific
  
  
  shared$type <- "both"
  pgc$type <- "PGC-specific"
  soma$type <- "soma-specific"
  
  all_boundaries <- append(append(shared, pgc), soma)
  all_boundaries$final_classification <- all_boundaries$type
  
  all_boundaries$final_classification[all_boundaries$type == "both" & (all_boundaries$score_in_PGC - all_boundaries$score_in_soma) > 0] <- "shared_stronger_in_PGC"
  all_boundaries$final_classification[all_boundaries$type == "both" & (all_boundaries$score_in_PGC - all_boundaries$score_in_soma) < 0] <- "shared_stronger_in_soma"
  
  all_boundaries_df <- as.data.frame(all_boundaries)
  
  assign(paste0("all_boundaries_df_", sex), all_boundaries_df[, c(1:3, 6:9)])
 
  p <- ggplot(data = all_boundaries_df) +
    geom_point(mapping = aes(x = score_in_soma, y = score_in_PGC, color = final_classification)) +
    theme_classic() + 
    ggtitle(paste("Boundaries in", sex)) +
    xlab("Boundary score in soma") + ylab("Boundary score in PGC") +
    scale_x_continuous(limits = c(0,2.5)) +   scale_y_continuous(limits = c(0,2.5)) + 
    geom_abline(linetype=2)  + coord_fixed(ratio = 1)
  
  
  ggsave(paste0("output/scatterplot_boundary_score_by_type_", sex , ".pdf"), useDingbats = FALSE,
         device = "pdf", width = 6, height = 5, p)  

  table_types_of_boundaries <- rbind(table_types_of_boundaries,
                                     table(all_boundaries_df$final_classification))
  rownames(table_types_of_boundaries)[nrow(table_types_of_boundaries)] <- sex
  
  

}

# PLOT table_types_of_boundaries
pdf(paste0("output/num_boundaries_by_type_", sex , ".pdf"), width = 6, height = 6)
par(mar = c(12,5,3,3))
bp <- barplot(table_types_of_boundaries, beside = T, col = c("dodgerblue4", "darkorange"),
              border = F,
              las = 2,
              ylab = "Number of boundaries",
              ylim = c(0, 4200),
              legend = T, args.legend = list(fill = c("darkorange", "dodgerblue4"),
                                             x = "topleft",
                                             legend = c("female", "male"),
                                             bty = "n", border = F))


text(labels = table_types_of_boundaries, 
     x = bp, y = table_types_of_boundaries, pos = 3, cex = 0.7)
dev.off()



###########


# Classification

for(sex in c("female", "male")) {
  all_boundaries <- get(paste0("all_boundaries_df_", sex))
  all_boundaries$percent_loss_inPGC <- (100- all_boundaries$score_in_PGC/all_boundaries$score_in_soma * 100)
  all_boundaries$log2_PGC_soma <- log2(all_boundaries$score_in_PGC/all_boundaries$score_in_soma)
  all_boundaries$log2_soma_PGC <- log2(all_boundaries$score_in_soma/all_boundaries$score_in_PGC)
  
  # Sort by change in strength
  all_boundaries <- all_boundaries[order(all_boundaries$log2_PGC_soma, decreasing = F),]
  
  all_boundaries$sex <- sex
  
  soma_stronger <- all_boundaries[all_boundaries$final_classification == "shared_stronger_in_soma", ]
  
  # I split them into 2 groups: 
  # 1) mild (strength in PGC is more than the half of score in soma)
  # 2) severe (strength in PGC is the half or less than the half of score in soma)
  soma_stronger_mild   <- soma_stronger[soma_stronger$log2_PGC_soma > -1, ]
  soma_stronger_severe <- soma_stronger[soma_stronger$log2_PGC_soma <= -1, ]
  soma_stronger_mild$final_classification <- "Mild_loss_in_PGC"
  soma_stronger_severe$final_classification <- "Severe_loss_in_PGC"
  
  # Other groups
  soma_specific <- all_boundaries[all_boundaries$final_classification == "soma-specific", ]
  soma_specific$final_classification <- "Total_loss_in_PGC"
  PGC_stronger <- all_boundaries[all_boundaries$final_classification == "shared_stronger_in_PGC", ]
  PGC_specific <- all_boundaries[all_boundaries$final_classification == "PGC-specific", ]
  
  # For the soma-specific and PGC-specific, since they all have the same score, I will randomise them. Otherwise,
  # due to their position sorted by genome, in some cases, some close boundaries create patterns in the heatmaps in which they enhance signal (or lack of)
  # because of the similarity between some lines that are have overlapping coordinates. This is corrected when putting them in random order. 
  set.seed(123)
  PGC_specific <- PGC_specific[sample(1:(nrow(PGC_specific)), size = nrow(PGC_specific), replace = F) , ]
  set.seed(123)
  soma_specific <- soma_specific[sample(1:(nrow(soma_specific)), size = nrow(soma_specific), replace = F) , ]

  
  # All boundaries in one table
  
  all <- rbind(PGC_specific, PGC_stronger, soma_stronger_mild, soma_stronger_severe, soma_specific)
  all <- all[, c(1,2,3,4,5,6,7)]
  write.table(all,
              col.names = F, quote = F, row.names = F, sep = "\t",
              file = paste0("output/classified_boundaries_", sex, "_mm39.bed"))
                            # Used to be called "all_boundaries_in_order_", sex, "_mm39.bed"
  
}

