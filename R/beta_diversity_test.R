# PERMANOVA for beta diversity
# test for beta-diversity (between groups diversity)
# global
# get data from phyloseq obj
# metadata <- as(sample_data(sasan_otu_final), "data.frame")
# metadata <- as(sample_data(full_OTU), "data.frame")
metadata <- data.frame(tidy_otu$samples)

counts <- as.matrix(tidy_otu$otu_table)
d <- vegan::vegdist(x = t(counts), method = "bray") # <-- distance between samples!
# d <- vegan::vegdist(x = counts, method = "bray")

permanova_stats <- adonis(formula = t(counts) ~ condition,
                          data = metadata,
                          method = "bray",
                          permutations = 1000)

# per day
# day 1
x1 <- metadata[which(metadata$cond_day == "aerobic1" | metadata$cond_day == "anaerobic1"), ]
y1 <- counts[, which(colnames(counts) %in% x1$filename)]
permanova_day1 <- adonis(formula = t(y1) ~ condition,
                         data = x1,
                         method = "bray",
                         permutations = 1000)

x1 <- metadata[which(metadata$cond_day == "aerobic2" | metadata$cond_day == "anaerobic2"), ]
y1 <- counts[, which(colnames(counts) %in% x1$filename)]
permanova_day2 <- adonis(formula = t(y1) ~ condition,
                         data = x1,
                         method = "bray",
                         permutations = 1000)

x1 <- metadata[which(metadata$cond_day == "aerobic3" | metadata$cond_day == "anaerobic3"), ]
y1 <- counts[, which(colnames(counts) %in% x1$filename)]
permanova_day3 <- adonis(formula = t(y1) ~ condition,
                         data = x1,
                         method = "bray",
                         permutations = 1000)


# PERMANOVA on human samples vs sasan's data
# tmp1 <- tidy_full$taxonomy$Genus
# tmp2 <- tidy_human$taxonomy$genus
# common_otus <- intersect(tmp1,tmp2)
# 
# id1 <- which(tidy_full$taxonomy$Genus %in% common_otus)
# id2 <- which(tidy_human$taxonomy$genus %in% common_otus)
# comb1 <- as.matrix(tidy_full$otu_table)[id1, ]
# comb2 <- as.matrix(tidy_human$otu_table)[id2, ]
# d1 <- as.vector(comb1)
# d2 <- as.vector(comb2)
# 
# otu_mat <- data_frame(genus = c(rep(tidy_full$taxonomy$Genus[id1], ncol(comb1)), 
#                                 rep(tidy_human$taxonomy$genus[id2], ncol(comb2))),
#                       counts = c(d1, d2),
#                       relative_abundance = c(as.vector(apply(comb1, 2, function(x) x / sum(x))),
#                                              as.vector(apply(comb2, 2, function(x) x / sum(x)))),
#                       minmax_counts = c(as.vector(apply(comb1, 2, function(x) (x - min(x)) / (max(x) - min(x)))),
#                                         as.vector(apply(comb2, 2, function(x) (x - min(x)) / (max(x) - min(x))))),
#                       scaled_data = c(as.vector(scale(comb1)),
#                                    as.vector(scale(comb2))),
#                       condition = c(as.vector(sapply(tidy_full$samples$condition, rep, nrow(comb1))),
#                                     rep("human_stool", length(d2))))
# 
# otu_mat %>% 
#   dplyr::filter(., counts != 0) %>% 
#   ggplot() + 
#   # geom_violin(aes(x = genus, y = log10(counts), fill = condition)) +
#   # geom_boxplot(aes(x = genus, y = log10(counts), fill = condition), outlier.alpha = 0) +
#   # geom_boxplot(aes(x = genus, y = log10(counts), fill = condition), notch = TRUE) +
#   # geom_boxplot(aes(x = genus, y = minmax_counts, fill = condition), outlier.alpha = 0) +
#   # geom_boxplot(aes(x = genus, y = scaled_data, fill = condition)) + 
#   # geom_point(aes(x = genus, y = log10(counts), color = condition), position = position_jitterdodge(jitter.width = 0.2)) +
#   geom_point(aes(x = genus, y = relative_abundance, color = condition), position = position_jitterdodge(jitter.width = 0.1), size = 3) +
#   # ylim(c(0, 0.5)) + 
#   # facet_grid(. ~ genus, scales = "free") +
#   # scale_fill_manual(values = c("#006699", "#00CCFF", "#FF3300"), 
#   #                   breaks = c("aerobic", "anaerobic", "human_stool"), 
#   #                   labels = c("Aerobic", "Anaerobic", "Human stool"), 
#   #                   name = "Sample type") + 
#   scale_color_manual(values = c("#006699", "#00CCFF", "#FF3300"), 
#                     breaks = c("aerobic", "anaerobic", "human_stool"), 
#                     labels = c("Aerobic", "Anaerobic", "Human stool"), 
#                     name = "Sample type") + 
#   labs(x = "Genus", y = "log10(counts)") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
#         axis.text.y = element_text(size = 12, color = "black"),
#         axis.title = element_text(size = 20, face = "bold", color = "black"),
#         panel.grid = element_blank())
#   
# 
