# richness analysis ----
message("Richness analysis...")

richness_data <- phyloseq::estimate_richness(full_OTU)
richness_data <- data_frame(metric = as.vector(sapply(colnames(richness_data), rep, nrow(richness_data))),
                            alpha_diversity = as.vector(as.matrix(richness_data)),
                            sample_name = as.vector(t(sapply(rownames(richness_data), rep, ncol(richness_data)))),
                            condition = as.vector(t(sapply(sample_data(full_OTU)[, 6], rep, ncol(richness_data)))),
                            replicate = as.vector(t(sapply(sample_data(full_OTU)[, 7], rep, ncol(richness_data)))),
                            day = as.vector(t(sapply(sample_data(full_OTU)[, 5], rep, ncol(richness_data)))),
                            cond_day = as.vector(t(sapply(sample_data(full_OTU)[, 8], rep, ncol(richness_data)))))

# add human data
rich_human <- phyloseq::estimate_richness(HMPv35)
stid <- sample_data(HMPv35)[which(sample_data(HMPv35)[, 6] == "Stool"), 1]
rich_human <- rich_human[which(gsub("X", "", rownames(rich_human)) %in% stid$X.SampleID), ]

richness_data <- data_frame(metric = as.vector(sapply(colnames(richness_data), rep, nrow(richness_data))),
                            alpha_diversity = as.vector(as.matrix(richness_data)),
                            sample_name = as.vector(t(sapply(rownames(richness_data), rep, ncol(richness_data)))),
                            condition = as.vector(t(sapply(sample_data(full_OTU)[, 6], rep, ncol(richness_data)))),
                            replicate = as.vector(t(sapply(sample_data(full_OTU)[, 7], rep, ncol(richness_data)))),
                            day = as.vector(t(sapply(sample_data(full_OTU)[, 5], rep, ncol(richness_data)))),
                            cond_day = as.vector(t(sapply(sample_data(full_OTU)[, 8], rep, ncol(richness_data)))))

combo_richness <- data_frame(metric = c(richness_data$metric, 
                                        as.vector(sapply(colnames(rich_human), rep, nrow(rich_human)))),
                             alpha_diversity = c(richness_data$alpha_diversity,
                                                 as.vector(as.matrix(rich_human))),
                             sample_name = c(richness_data$sample_name, 
                                             rep("human_stool", nrow(rich_human) * ncol(rich_human))),
                             condition = c(richness_data$condition, 
                                           rep("human_stool", nrow(rich_human) * ncol(rich_human))),
                             replicate = c(richness_data$replicate,
                                           rep(1, nrow(rich_human) * ncol(rich_human))),
                             day = c(richness_data$day,
                                     rep(NA, nrow(rich_human) * ncol(rich_human))),
                             cond_day = c(richness_data$cond_day,
                                          rep(1, nrow(rich_human) * ncol(rich_human))))




