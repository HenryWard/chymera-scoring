#!/usr/bin/env Rscript

library(ggplot2)
library(ggthemes)

######
# PARAMETER SETTING
######

# MODIFY THE THREE VARIABLES BELOW TO RUN THIS SCRIPT ON YOUR MACHINE!!!

# Sets working directory, input file location, output folder path, and cell line.
# Guide type must be one of "single", "dual", "paralog_single" or "paralog_dual"
setwd("/project/csbio/henry/Documents/projects/dual_guide/chymera_github")
output_folder <- file.path("output", "dual_targeted")
guide_type <- "dual"

######
# MAIN SCRIPT
######

# Makes output folder if it doesn't exist
if(!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Optionally removes gene pairs with too few/many guid
rm_guides <- TRUE

# Loads in data
input_file_paralogs <- file.path("input", "torin_norm_counts.txt")
input_file_dual_targeted <- file.path("input", "dual_targeted_counts.txt")
paralogs <- read.csv(input_file_paralogs, sep = "\t", stringsAsFactors = FALSE)
torin <- read.csv(input_file_dual_targeted, sep = "\t", stringsAsFactors = FALSE)

# Calculates LFCs for both timepoints and cell lines for paralog data
paralogs$HAP1_T0 <- log2(paralogs$HAP1_T0 + 0.1)
paralogs$HAP1_T12_logFCA <- log2(paralogs$HAP1_T12A + 0.1) - paralogs$HAP1_T0 
paralogs$HAP1_T12_logFCB <- log2(paralogs$HAP1_T12B + 0.1) - paralogs$HAP1_T0 
paralogs$HAP1_T12_logFCC <- log2(paralogs$HAP1_T12C + 0.1) - paralogs$HAP1_T0 
paralogs$HAP1_T18_logFCA <- log2(paralogs$HAP1_T18A + 0.1) - paralogs$HAP1_T0 
paralogs$HAP1_T18_logFCB <- log2(paralogs$HAP1_T18B + 0.1) - paralogs$HAP1_T0 
paralogs$HAP1_T18_logFCC <- log2(paralogs$HAP1_T18C + 0.1) - paralogs$HAP1_T0
paralogs$HAP1_Torin_T12_logFCA <- log2(paralogs$HAP1_Torin_T12A + 0.1) - paralogs$HAP1_T0 
paralogs$HAP1_Torin_T12_logFCB <- log2(paralogs$HAP1_Torin_T12B + 0.1) - paralogs$HAP1_T0 
paralogs$HAP1_Torin_T12_logFCC <- log2(paralogs$HAP1_Torin_T12C + 0.1) - paralogs$HAP1_T0 
paralogs$HAP1_Torin_T18_logFCA <- log2(paralogs$HAP1_Torin_T18A + 0.1) - paralogs$HAP1_T0 
paralogs$HAP1_Torin_T18_logFCB <- log2(paralogs$HAP1_Torin_T18B + 0.1) - paralogs$HAP1_T0 
paralogs$HAP1_Torin_T18_logFCC <- log2(paralogs$HAP1_Torin_T18C + 0.1) - paralogs$HAP1_T0

# Parses out paralog single and dual-targeting pairs
paralog_singles <- paralogs[paralogs$Gene_symbol1 == "---" | paralogs$Gene_symbol2 == "---",]
paralog_singles <- paralog_singles[paralog_singles$Library != "DualTargeting_iCtrl",]
paralog_singles <- paralog_singles[paralog_singles$Cas9_Guide_Type != paralog_singles$Cpf1_Guide_Type, ]
paralog_singles$Gene_symbol1[paralog_singles$Gene_symbol1 == "---"] <- paralog_singles$Gene_symbol2[paralog_singles$Gene_symbol1 == "---"]
paralog_dual <- paralogs[paralogs$Library == "Paralogs",]
paralog_dual <- paralog_dual[paralog_dual$Gene_symbol1 != "---" & paralog_dual$Gene_symbol2 != "---",]

# Gets all unique paralog pairs
gene_pairs <- t(unlist(apply(paralog_dual[,c("Gene_symbol1", "Gene_symbol2")], 1, function(x) sort(x))))
gene_pairs <- gene_pairs[!duplicated(gene_pairs),]
gene_pairs <- gene_pairs[gene_pairs[,1] %in% paralog_singles$Gene_symbol1 & gene_pairs[,2] %in% paralog_singles$Gene_symbol1,]
unique_paralog_genes <- unique(c(gene_pairs[,1], gene_pairs[,2]))
paralog_singles <- paralog_singles[paralog_singles$Gene_symbol1 %in% unique_paralog_genes,]

# Calculates LFCs for both timepoints and cell lines for non-paralog data
colnames(torin)[colnames(torin) == "Gene.symbol1"] <- "Gene_symbol1"
torin$HAP1.T0 <- log2(torin$HAP1.T0 + 0.1)
torin$HAP1_T12_logFC <- rowMeans(torin[,c("HAP1.T12A", "HAP1.T12B", "HAP1.T12C")], na.rm = TRUE)
torin$HAP1_T12_logFC <- log2(torin$HAP1_T12_logFC + 0.1) -  torin$HAP1.T0
torin$HAP1_T12_logFCA <- log2(torin$HAP1.T12A + 0.1) -  torin$HAP1.T0 
torin$HAP1_T12_logFCB <- log2(torin$HAP1.T12B + 0.1) -  torin$HAP1.T0
torin$HAP1_T12_logFCC <- log2(torin$HAP1.T12C + 0.1) -  torin$HAP1.T0 
torin$HAP1_T18_logFC <- rowMeans(torin[,c("HAP1.T18A", "HAP1.T18B", "HAP1.T18C")], na.rm = TRUE)
torin$HAP1_T18_logFC <- log2(torin$HAP1_T18_logFC + 0.1) -  torin$HAP1.T0 
torin$HAP1_T18_logFCA <- log2(torin$HAP1.T18A + 0.1) -  torin$HAP1.T0 
torin$HAP1_T18_logFCB <- log2(torin$HAP1.T18B + 0.1) -  torin$HAP1.T0
torin$HAP1_T18_logFCC <- log2(torin$HAP1.T18C + 0.1) -  torin$HAP1.T0 
torin$HAP1_Torin_T12_logFC <- rowMeans(torin[,c("HAP1.Torin.T12A", "HAP1.Torin.T12B", "HAP1.Torin.T12C")], na.rm = TRUE)
torin$HAP1_Torin_T12_logFC <- log2(torin$HAP1_Torin_T12_logFC + 0.1) -  torin$HAP1.T0
torin$HAP1_Torin_T12_logFCA <- log2(torin$HAP1.Torin.T12A + 0.1) -  torin$HAP1.T0 
torin$HAP1_Torin_T12_logFCB <- log2(torin$HAP1.Torin.T12B + 0.1) -  torin$HAP1.T0 
torin$HAP1_Torin_T12_logFCC <- log2(torin$HAP1.Torin.T12C + 0.1) -  torin$HAP1.T0 
torin$HAP1_Torin_T18_logFC <- rowMeans(torin[,c("HAP1.Torin.T18A", "HAP1.Torin.T18B", "HAP1.Torin.T18C")], na.rm = TRUE)
torin$HAP1_Torin_T18_logFC <- log2(torin$HAP1_Torin_T18_logFC + 0.1) -  torin$HAP1.T0
torin$HAP1_Torin_T18_logFCA <- log2(torin$HAP1.Torin.T18A + 0.1) -  torin$HAP1.T0  
torin$HAP1_Torin_T18_logFCB <- log2(torin$HAP1.Torin.T18B + 0.1) -  torin$HAP1.T0  
torin$HAP1_Torin_T18_logFCC <- log2(torin$HAP1.Torin.T18C + 0.1) -  torin$HAP1.T0  

# Separates out dual-targeted guides from single-targeted guides
dual <- torin[torin$Library == "DualTargeting",]
dual <- dual[dual$Cas12a != "inter" & dual$Cas9 != "inter",]
dual <- dual[dual$Gene_symbol1 != "---" & dual$Gene_symbol1 != "NT",]
singles <- torin[torin$Cas9 == "inter" | torin$Cas12a == "inter",]
singles <- singles[singles$Cas9 != singles$Cas12a,]
singles <- singles[singles$Gene_symbol1 != "---" & singles$Gene_symbol1 != "NT",]

# Gets residual LFCs
singles$early_residual_logFC <- singles$HAP1_Torin_T12_logFC - singles$HAP1_T12_logFC
singles$late_residual_logFC <- singles$HAP1_Torin_T18_logFC - singles$HAP1_T18_logFC
dual$early_residual_logFC <- dual$HAP1_Torin_T12_logFC - dual$HAP1_T12_logFC
dual$late_residual_logFC <- dual$HAP1_Torin_T18_logFC - dual$HAP1_T18_logFC
paralog_singles$early_residual_logFC <-
  paralog_singles$HAP1_Torin_T12_logFC - paralog_singles$HAP1_T12_logFC
paralog_singles$late_residual_logFC <-
  paralog_singles$HAP1_Torin_T18_logFC - paralog_singles$HAP1_T18_logFC

# Calculates log FCs and p-values for paralog single genes
df <- NA
if (guide_type == "single") {
  df <- data.frame(gene = unique(singles$Gene_symbol1), stringsAsFactors = FALSE)
} else if (guide_type == "dual") {
  df <- data.frame(gene = unique(dual$Gene_symbol1), stringsAsFactors = FALSE)
} else if (guide_type == "paralog_single") {
  df <- data.frame(gene = unique_paralog_genes, stringsAsFactors = FALSE)
} else if (guide_type == "paralog_dual") {
  df <- data.frame(gene = paste(gene_pairs[,1], gene_pairs[,2], sep = "/"), stringsAsFactors = FALSE)
}
stat_table <- data.frame(gene = df$gene, 
                         hap1_time1_n_observed = rep(NA, length(df$gene)), 
                         hap1_time2_n_observed = rep(NA, length(df$gene)),
                         torin_time1_n_observed = rep(NA, length(df$gene)), 
                         torin_time2_n_observed = rep(NA, length(df$gene)))
results_file <- file.path(output_folder, "results.txt")
df$hap1_time1_observed_logFC <- NA
df$torin_time1_observed_logFC <- NA
df$time1_pval <- NA
df$hap1_time2_observed_logFC <- NA
df$torin_time2_observed_logFC <- NA
df$time2_pval <- NA
hap1_time1_observed_sets <- list()
torin_time1_observed_sets <- list()
hap1_time2_observed_sets <- list()
torin_time2_observed_sets <- list()
for (i in 1:nrow(df)) {
  gene <- df$gene[i]
  torin_time1_observed <- c()
  torin_time2_observed <- c()
  hap1_time1_observed <- c()
  hap1_time2_observed <- c()
  if (guide_type == "single") {
    temp <- singles[singles$Gene_symbol1 == gene,]
    torin_time1_observed <- c(temp$HAP1_Torin_T12_logFCA, temp$HAP1_Torin_T12_logFCB, temp$HAP1_Torin_T12_logFCC)
    torin_time2_observed <- c(temp$HAP1_Torin_T18_logFCA, temp$HAP1_Torin_T18_logFCB, temp$HAP1_Torin_T18_logFCC)
    hap1_time1_observed <- c(temp$HAP1_T12_logFCA, temp$HAP1_T12_logFCB, temp$HAP1_T12_logFCC)
    hap1_time2_observed <- c(temp$HAP1_T18_logFCA, temp$HAP1_T18_logFCB, temp$HAP1_T18_logFCC)
    if (rm_guides) {
      if (length(torin_time1_observed) != 15) { next }
    }
  } else if (guide_type == "dual") {
    temp <- dual[dual$Gene_symbol1 == gene,]
    torin_time1_observed <- c(temp$HAP1_Torin_T12_logFCA, temp$HAP1_Torin_T12_logFCB, temp$HAP1_Torin_T12_logFCC)
    torin_time2_observed <- c(temp$HAP1_Torin_T18_logFCA, temp$HAP1_Torin_T18_logFCB, temp$HAP1_Torin_T18_logFCC)
    hap1_time1_observed <- c(temp$HAP1_T12_logFCA, temp$HAP1_T12_logFCB, temp$HAP1_T12_logFCC)
    hap1_time2_observed <- c(temp$HAP1_T18_logFCA, temp$HAP1_T18_logFCB, temp$HAP1_T18_logFCC)
    if (rm_guides) {
      if (length(torin_time1_observed) != 18) { 
        next 
      } else { 
        torin_time1_observed <- sample(torin_time1_observed, 15)
        torin_time2_observed <- sample(torin_time2_observed, 15)
        hap1_time1_observed <- sample(hap1_time1_observed, 15)
        hap1_time2_observed <- sample(hap1_time2_observed, 15)
      }
    }
  } else if (guide_type == "paralog_single") {
    temp <- paralog_singles[paralog_singles$Gene_symbol1 == gene,]
    torin_time1_observed <- c(temp$HAP1_Torin_T12_logFCA, temp$HAP1_Torin_T12_logFCB, temp$HAP1_Torin_T12_logFCC)
    torin_time2_observed <- c(temp$HAP1_Torin_T18_logFCA, temp$HAP1_Torin_T18_logFCB, temp$HAP1_Torin_T18_logFCC)
    hap1_time1_observed <- c(temp$HAP1_T12_logFCA, temp$HAP1_T12_logFCB, temp$HAP1_T12_logFCC)
    hap1_time2_observed <- c(temp$HAP1_T18_logFCA, temp$HAP1_T18_logFCB, temp$HAP1_T18_logFCC)
    if (rm_guides) {
      if (length(torin_time1_observed) != 24) { next }
    }
  } else if (guide_type == "paralog_dual") {
    ind1 <- paralog_dual$Gene_symbol1 == gene_pairs[i,1] & paralog_dual$Gene_symbol2 == gene_pairs[i,2]
    ind2 <- paralog_dual$Gene_symbol1 == gene_pairs[i,2] & paralog_dual$Gene_symbol2 == gene_pairs[i,1]
    temp <- paralog_dual[ind1 | ind2,]
    if (rm_guides) {
      if (nrow(temp) < 24) { 
        next 
      } else {
        temp <- temp[sample(rownames(temp), 24),]
      }
    }
    torin_time1_observed <- temp$HAP1_Torin_T12_logFC
    torin_time2_observed <- temp$HAP1_Torin_T18_logFC
    hap1_time1_observed <- temp$HAP1_T12_logFC
    hap1_time2_observed <- temp$HAP1_T18_logFC
  }
  stat_table$torin_time1_n_observed[i] <- length(torin_time1_observed)
  stat_table$torin_time2_n_observed[i] <- length(torin_time2_observed)
  stat_table$hap1_time1_n_observed[i] <- length(hap1_time1_observed)
  stat_table$hap1_time2_n_observed[i] <- length(hap1_time2_observed)
  df$hap1_time1_observed_logFC[i] <- mean(c(hap1_time1_observed), na.rm = TRUE)
  df$hap1_time2_observed_logFC[i] <- mean(c(hap1_time2_observed), na.rm = TRUE)
  df$torin_time1_observed_logFC[i] <- mean(c(torin_time1_observed), na.rm = TRUE)
  df$torin_time2_observed_logFC[i] <- mean(c(torin_time2_observed), na.rm = TRUE)
  df$time1_pval[i] <- suppressWarnings(wilcox.test(hap1_time1_observed, torin_time1_observed)$p.value)
  df$time2_pval[i] <- suppressWarnings(wilcox.test(hap1_time2_observed, torin_time2_observed)$p.value)
  hap1_time1_observed_sets[[i]] <- hap1_time1_observed
  torin_time1_observed_sets[[i]] <- torin_time1_observed
  hap1_time2_observed_sets[[i]] <- hap1_time2_observed
  torin_time2_observed_sets[[i]] <- torin_time2_observed
}
stat_table <- stat_table[complete.cases(stat_table),]
df <- df[complete.cases(df),]
df$time1_residual_logFC <- df$torin_time1_observed_logFC - df$hap1_time1_observed_logFC
df$time2_residual_logFC <- df$torin_time2_observed_logFC - df$hap1_time2_observed_logFC
df <- df[!is.na(df$time1_residual_logFC) & !is.na(df$time2_residual_logFC),]

# Writes stat table to file
write.table(stat_table, file.path(output_folder, "stat_table.tsv"), sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)

# Calls significant differences with regards to torin response
df$time1_fdr <- p.adjust(df$time1_pval, "BH", nrow(df))
df$time2_fdr <- p.adjust(df$time2_pval, "BH", nrow(df))
df$time1_significant <- df$time1_fdr < 0.01
df$time2_significant <- df$time2_fdr < 0.01

# Calls suppressors and sensitizers
df$time1_response <- "None"
df$time1_response[df$time1_significant & df$time1_residual_logFC < 0] <- "Sensitizer"
df$time1_response[df$time1_significant & df$time1_residual_logFC > 0] <- "Suppressor"
df$time2_response <- "None"
df$time2_response[df$time2_significant & df$time2_residual_logFC < 0] <- "Sensitizer"
df$time2_response[df$time2_significant & df$time2_residual_logFC > 0] <- "Suppressor"

# Prints number of tentative significant guides
sink(file.path(output_folder, "stats.txt"))
cat(paste0("Number of early genes with p-val < 0.01: ", sum(df$time1_pval < 0.01), "\n"))
cat(paste0("Number of late genes with p-val < 0.01: ", sum(df$time2_pval < 0.01), "\n"))
cat(paste0("Number of early genes with FDR < 0.01: ", sum(df$time1_significant), "\n"))
cat(paste0("Number of late genes with FDR < 0.01: ", sum(df$time2_significant), "\n"))
sink()

# Plots HAP1 with and without torin
if (sum(df$time1_significant) > 0) {
  ggplot(df, aes(x = hap1_time1_observed_logFC, y = torin_time1_observed_logFC)) +
    geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    geom_vline(xintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    geom_abline(slope = 1, intercept = 0, size = 1.5, alpha = 0.5, color = "Black") +
    geom_point(aes(color = time1_response)) +
    scale_color_manual(values = c("Gray", "Blue", "Yellow")) +
    xlab("HAP1 - torin T18 mean log FC") +
    ylab("HAP1 + torin T18 mean log FC") +
    labs(color = "Significant torin response") +
    theme_tufte(base_size = 14) +
    theme(legend.position = c(0.15, 0.9),
          axis.text.x = element_text(color = "Black", size = 14),
          axis.text.y = element_text(color = "Black", size = 14),
          legend.text = element_text(size = 14))
  ggsave(file.path(output_folder, "early_differential.png"))
}
if (sum(df$time2_significant) > 0) {
  ggplot(df, aes(x = hap1_time2_observed_logFC, y = torin_time2_observed_logFC)) +
    geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    geom_vline(xintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    geom_abline(slope = 1, intercept = 0, size = 1.5, alpha = 0.5, color = "Black") +
    geom_point(aes(color = time2_response)) +
    scale_color_manual(values = c("Gray", "Blue", "Yellow")) +
    xlab("HAP1 - torin T18 mean log FC") +
    ylab("HAP1 + torin T18 mean log FC") +
    labs(color = "Significant torin response") +
    theme_tufte(base_size = 14) +
    theme(legend.position = c(0.15, 0.9),
          axis.text.x = element_text(color = "Black", size = 14),
          axis.text.y = element_text(color = "Black", size = 14),
          legend.text = element_text(size = 14))
  ggsave(file.path(output_folder, "late_differential.png"))
}

# Makes rank plots
n_labels <- 5
temp_df <- cbind(df[order(df$time1_residual_logFC),], data.frame(ind = 1:nrow(df)))
temp_df1 <- temp_df[temp_df$time1_response == "None",]
temp_df2 <- temp_df[temp_df$time1_response != "None",]
early_l_labels <- temp_df1[order(temp_df1$time1_residual_logFC, decreasing = FALSE),]
early_l_labels <- paste(rev(early_l_labels[1:n_labels,"gene"]), collapse = "\n")
early_r_labels <- temp_df1[order(temp_df1$time1_residual_logFC, decreasing = TRUE),]
early_r_labels <- paste(early_r_labels[1:n_labels,"gene"], collapse = "\n")
late_l_labels <- temp_df1[order(temp_df1$time2_residual_logFC, decreasing = FALSE),]
late_l_labels <- paste(rev(late_l_labels[1:n_labels,"gene"]), collapse = "\n")
late_r_labels <- temp_df1[order(temp_df1$time2_residual_logFC, decreasing = TRUE),]
late_r_labels <- paste(late_r_labels[1:n_labels,"gene"], collapse = "\n")
ggplot() +
  geom_hline(yintercept = 0, linetype = 2, size = 0.8, alpha = 1, color = "Gray") +
  geom_point(data = temp_df1, aes(x = ind, y = time1_residual_logFC, color = time1_response)) +
  geom_point(data = temp_df2, aes(x = ind, y = time1_residual_logFC, fill = time1_response),
             size = 2.5, pch = 21) +
  #annotate("text", x = 400, y = -1.7, label = early_l_labels) +
  #annotate("text", x = 4200, y = 2.82, label = early_r_labels) +
  scale_color_manual(values = c("Gray")) +
  scale_fill_manual(values = c("Blue", "Yellow")) +
  xlab("Rank") +
  ylab("Mean residual log FC") +
  labs(fill = "Significant torin response") +
  guides(color = FALSE) +
  theme_tufte(base_size = 16) +
  theme(legend.position = c(0.25, 0.9),
        axis.text.x = element_text(color = "Black", size = 14),
        axis.text.y = element_text(color = "Black", size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18))
ggsave(file.path(output_folder, "early_gene_rank.png"))
temp_df <- cbind(df[order(df$time2_residual_logFC),], data.frame(ind = 1:nrow(df)))
temp_df1 <- temp_df[temp_df$time2_response == "None",]
temp_df2 <- temp_df[temp_df$time2_response != "None",]
ggplot() +
  geom_hline(yintercept = 0, linetype = 2, size = 0.8, alpha = 1, color = "Gray") +
  geom_point(data = temp_df1, aes(x = ind, y = time2_residual_logFC, color = time2_response)) +
  geom_point(data = temp_df2, aes(x = ind, y = time2_residual_logFC, fill = time2_response),
             size = 2.5, pch = 21) +
  #annotate("text", x = 450, y = -2.5, label = late_l_labels) +
  #annotate("text", x = 4250, y = 4.88, label = late_r_labels) +
  scale_color_manual(values = c("Gray")) +
  scale_fill_manual(values = c("Blue", "Yellow")) +
  xlab("Rank") +
  ylab("Mean residual log FC") +
  labs(fill = "Significant torin response") +
  guides(color = FALSE) +
  theme_tufte(base_size = 16) +
  theme(legend.position = c(0.25, 0.9),
        axis.text.x = element_text(color = "Black", size = 14),
        axis.text.y = element_text(color = "Black", size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18))
ggsave(file.path(output_folder, "late_gene_rank.png"))

# Writes sorted data to file
df_time1 <- df[,c(1, which(grepl("time1", colnames(df))))]
df_time2 <- df[,c(1, which(grepl("time2", colnames(df))))] 
df_time1 <- df_time1[order(df_time1$time1_residual_logFC, decreasing = TRUE),]
df_time2 <- df_time2[order(df_time2$time2_residual_logFC, decreasing = TRUE),]
write.table(df, file.path(output_folder, "differential_calls.tsv"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(df_time1, file.path(output_folder, "early_calls.tsv"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(df_time2, file.path(output_folder, "late_calls.tsv"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Writes lists of guides to file
saveRDS(hap1_time1_observed_sets, file.path(output_folder, "early_hap1.rds"))
saveRDS(hap1_time2_observed_sets, file.path(output_folder, "late_hap1.rds"))
saveRDS(torin_time1_observed_sets, file.path(output_folder, "early_torin.rds"))
saveRDS(torin_time2_observed_sets, file.path(output_folder, "late_torin.rds"))
