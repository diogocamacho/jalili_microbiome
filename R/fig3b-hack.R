X <- import_biom("/Volumes/THoR/omics/data/microbial_community_sequencing/organ_chip/processed/media_optimization/otu_table-json.biom")


# get alessio's samples
f24_samples <- grep("F24",colnames(otu_table(X)))
hs_stock <- grep("Hmb1$",colnames(otu_table(X)))
aOTU <- prune_samples(samples = colnames(otu_table(X))[c(f24_samples,hs_stock)],X)

cols <- gsub(".F24","",gsub("Wyss.","",colnames(otu_table(aOTU))))
sample_id <- c(1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,"Hmb1")
media <- c(rep("20% FBS 4g/L glucose",3),rep("Mal-suco 0.12 g/L",3),rep("Pectin 0.5 g/L",3),rep("Starch 0.5 g/L",3),rep("Mucus 0.5 g/l",3),"Mal-suco 0.12 g/L","Human microbiome stock")

sample_metadata <- data.frame(sample_names=colnames(otu_table(aOTU)),
                              sample_id=cols,
                              media=media[match(cols,sample_id)])
rownames(sample_metadata) <- sample_metadata$sample_names
sample_data(aOTU) <- sample_metadata

xx <- tidy_phyloseq(phyloseq_obj = aOTU)

a1 <- rownames(tidy_full$otu_table)
a2 <- rownames(xx$otu_table)
xx2 <- xx$otu_table[match(a1, a2), ]
xx2 <- xx2[, which(xx$samples$sample_id == "Hmb1")]

tidy_full2 <- tidy_full
tidy_full2$otu_table <- tidy_full2$otu_table[, -25]
tidy_full2$samples <- tidy_full2$samples[-25, ]

a3 <- as_data_frame(tidy_full2$samples)
tidy_full2$samples <- as_data_frame(rbind(as.matrix(a3),
      c("Wyss.Hmb1", "Hmb1", "Human microbiome stock", NA, NA, NA, NA, NA)))


M <- tidy_correct_counts(tidy_phylo = tidy_full2,control_sample_id = 25)
summ_full <- tidy_taxa(tidy_phylo = M,summary_level = "genus")
relative_full <- relative_abundances(summary_matrix = summ_full$counts)

df4 <- data_frame(sample_name = as.vector(sapply(M$samples$filename, rep, nrow(summ_full$counts))),
                  sample_group = as.vector(sapply(M$samples$condition, rep, nrow(summ_full$counts))),
                  day=as.vector(sapply(M$samples$day,rep,nrow(summ_full$counts))),
                  replicate=as.vector(sapply(M$samples$replicate,rep,nrow(summ_full$counts))),
                  genus=rep(summ_full$bacteria,ncol(summ_full$counts)),
                  counts=as.vector(summ_full$counts),
                  relative_abundance=as.vector(relative_full))#,
# condition=as.vector(sapply(M$samples$condition,rep,nrow(summ_full$counts))),
# plot_var=paste(condition,"-- day",day,"-- replicate",replicate))
df4$replicate <- as.numeric(df4$replicate)

df4$sample_name[is.na(df4$sample_name)] <- "Hmb1"
df4$sample_group[is.na(df4$sample_group)] <- "Hmb1"
df4$day[is.na(df4$day)] <- "Hmb1"


df4 %>%
  # tibble::add_column(., s_id = factor(df2$sample_group, levels = c("Hmb11", "anaerobic", "aerobic"))) %>% 
  tibble::add_column(., p_id = gsub("Hmb1 \\w* \\(replicate NA\\)", "Hmb1", gsub(" )", ")", paste(df4$sample_group, df4$day, "(replicate", df4$replicate, ")")))) %>%
  ggplot() + 
  geom_bar(aes(x = p_id, y = relative_abundance, fill = genus), stat = "identity", color = "black") + 
  scale_fill_manual(breaks = c("Unknown", "Sutterella", "Bilophila", "Akkermansia", "Blautia",
                               "Oscillospira", "Ruminococcus", "Enterococcus", "[Eubacterium]",
                               "Parabacteroides", "Bacteroides", "Citrobacter"),
                    values = c("#CCCCCC", "#FFCC00", "#99CCCC", "#33CC33", "#CC3399",
                               "#CC99CC", "#336699", "#66CCFF", "#669933",
                               "#FF6600", "#993300", "#CCCCFF"),
                    name = "Genus") + 
  labs(x = "Sample", y = "Relative\nabundance") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 20, face = "bold", color = "black"),
        axis.text.x = element_text(size = 12, angle = 90, color = "black", hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12, color = "black"),
        panel.grid = element_blank())
ggsave(file=paste0("./results/Camacho_FigS6_",format(Sys.Date(),"%Y-%m-%d"),".pdf"), 
       width = 11, 
       height = 8, 
       units = "in", 
       dpi = 600)