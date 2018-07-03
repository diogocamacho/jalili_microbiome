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


