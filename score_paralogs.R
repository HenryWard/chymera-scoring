#!/usr/bin/env Rscript

library(ggplot2)
library(ggthemes)

######
# PARAMETER SETTING
######

# MODIFY THE FOUR VARIABLES BELOW TO RUN THIS SCRIPT ON YOUR MACHINE!!!

# Sets working directory, input file location, output folder path, and cell line.
# Cell line must be one of "hap1", "rpe1", or "hap1_torin"
setwd("/project/csbio/henry/Documents/projects/dual_guide/chymera_github")
input_file <- file.path("..", "input", "torin_norm_counts.txt")
output_folder <- file.path("output", "paralog_rpe1")
which_cell_line <- "rpe1"

######
# MAIN SCRIPT
######

# Makes output folder if it doesn't exist
if(!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Loads in data
torin <- read.csv(file.path(input_file), sep = "\t", stringsAsFactors = FALSE)
paralogs <- torin[torin$Library == "Paralogs",]
paralogs <- paralogs[,c("Gene_symbol1", "Gene_symbol2", 
                        "HAP1_Torin_T12_logFC", "HAP1_Torin_T18_logFC",
                        "HAP1_T12_logFC", "HAP1_T18_logFC",
                        "RPE1_T18_logFC", "RPE1_T24_logFC")]
singles <- torin[torin$Gene_symbol1 == "---" | torin$Gene_symbol2 == "---",]
singles <- singles[singles$Library != "DualTargeting_iCtrl",]
paralogs <- paralogs[!(rownames(paralogs) %in% rownames(singles)),]
gene_pairs <- paralogs[!duplicated(paralogs[,c("Gene_symbol1", "Gene_symbol2")]), c("Gene_symbol1", "Gene_symbol2")]
df <- data.frame(gene1 = gene_pairs$Gene_symbol1, gene2 = gene_pairs$Gene_symbol2, stringsAsFactors = FALSE)

# Calculates log FCs and p-values on an orientation level
stat_table <- data.frame(gene1 = df$gene1, gene2 = df$gene2, 
                         time1_n_observed = rep(NA, length(df$gene1)), time2_n_observed = rep(NA, length(df$gene1)),
                         n_expected_gene1_cas9 = rep(NA, length(df$gene1)), n_expected_gene2_cpf1 = rep(NA, length(df$gene1)))
results_file <- file.path(output_folder, paste0(which_cell_line, "_results.txt"))
sink(results_file)
df$time1_observed_logFC <- NA
df$time1_expected_logFC <- NA
df$time1_pval <- NA
df$time2_observed_logFC <- NA
df$time2_expected_logFC <- NA
df$time2_pval <- NA
time1_expected_sets <- list()
time1_observed_sets <- list()
time2_expected_sets <- list()
time2_observed_sets <- list()
for (i in 1:nrow(df)) {
  gene1 <- df$gene1[i]
  gene2 <- df$gene2[i]
  time1_all_expected <- c()
  time1_all_observed <- c()
  time2_all_expected <- c()
  time2_all_observed <- c()
  expected_guides <- singles[singles$Gene_symbol1 == gene1 | singles$Gene_symbol2 == gene2,] 
  stat_table$n_expected_gene1_cas9[i] <- max(nrow(expected_guides[expected_guides$Gene_symbol1 == gene1,]), 0, na.rm = TRUE)
  stat_table$n_expected_gene2_cpf1[i] <- max(nrow(expected_guides[expected_guides$Gene_symbol2 == gene2,]), 0, na.rm = TRUE)
  possible_pairs <- combn(rownames(expected_guides), 2)
  for (j in 1:ncol(possible_pairs)) {
    row1 <- singles[possible_pairs[1,j],]
    row2 <- singles[possible_pairs[2,j],]
    if (row1$Gene_symbol2 == row2$Gene_symbol1 | row1$Gene_symbol1 == row2$Gene_symbol2) {
      #print(rbind(row1, row2)[,c("Gene_symbol1", "Gene_symbol2", "HAP1_Torin_T12_logFC")])
      if (which_cell_line == "hap1_torin") {
        time1_all_expected <- c(time1_all_expected, row1$HAP1_Torin_T12_logFC + row2$HAP1_Torin_T12_logFC) 
        time2_all_expected <- c(time2_all_expected, row1$HAP1_Torin_T18_logFC + row2$HAP1_Torin_T18_logFC) 
      } else if (which_cell_line == "hap1") {
        time1_all_expected <- c(time1_all_expected, row1$HAP1_T12_logFC + row2$HAP1_T12_logFC) 
        time2_all_expected <- c(time2_all_expected, row1$HAP1_T18_logFC + row2$HAP1_T18_logFC) 
      } else if (which_cell_line == "rpe1") {
        time1_all_expected <- c(time1_all_expected, row1$RPE1_T18_logFC + row2$RPE1_T18_logFC) 
        time2_all_expected <- c(time2_all_expected, row1$RPE1_T24_logFC + row2$RPE1_T24_logFC) 
      }
    }
  }
  if (which_cell_line == "hap1_torin") {
    time1_all_observed <- paralogs$HAP1_Torin_T12_logFC[paralogs$Gene_symbol1 == gene1 & paralogs$Gene_symbol2 == gene2]
    time2_all_observed <- paralogs$HAP1_Torin_T18_logFC[paralogs$Gene_symbol1 == gene1 & paralogs$Gene_symbol2 == gene2]
  } else if (which_cell_line == "hap1") {
    time1_all_observed <- paralogs$HAP1_T12_logFC[paralogs$Gene_symbol1 == gene1 & paralogs$Gene_symbol2 == gene2]
    time2_all_observed <- paralogs$HAP1_T18_logFC[paralogs$Gene_symbol1 == gene1 & paralogs$Gene_symbol2 == gene2]
  } else if (which_cell_line == "rpe1") {
    time1_all_observed <- paralogs$RPE1_T18_logFC[paralogs$Gene_symbol1 == gene1 & paralogs$Gene_symbol2 == gene2]
    time2_all_observed <- paralogs$RPE1_T24_logFC[paralogs$Gene_symbol1 == gene1 & paralogs$Gene_symbol2 == gene2]
  }
  stat_table$time1_n_observed[i] <- length(time1_all_observed)
  stat_table$time2_n_observed[i] <- length(time2_all_observed)
  df$time1_observed_logFC[i] <- mean(c(time1_all_observed), na.rm = TRUE)
  df$time1_expected_logFC[i] <- mean(c(time1_all_expected), na.rm = TRUE)
  df$time2_observed_logFC[i] <- mean(c(time2_all_observed), na.rm = TRUE)
  df$time2_expected_logFC[i] <- mean(c(time2_all_expected), na.rm = TRUE)
  df$time1_pval[i] <- suppressWarnings(wilcox.test(time1_all_observed, time1_all_expected)$p.value)
  df$time2_pval[i] <- suppressWarnings(wilcox.test(time2_all_observed, time2_all_expected)$p.value)
  cat(paste(df$gene1[i], "and", df$gene2[i], "have", length(time1_all_expected), "expected and", 
            length(time1_all_observed), "observed\n"))
  time1_expected_sets[[i]] <- time1_all_expected
  time1_observed_sets[[i]] <- time1_all_observed
  time2_expected_sets[[i]] <- time2_all_expected
  time2_observed_sets[[i]] <- time2_all_observed
}
df <- df[complete.cases(df),]
df$time1_residual_logFC <- df$time1_observed_logFC - df$time1_expected_logFC
df$time2_residual_logFC <- df$time2_observed_logFC - df$time2_expected_logFC
df <- df[!is.na(df$time1_residual_logFC) & !is.na(df$time2_residual_logFC),]

# Writes stat table to file
write.table(stat_table, file.path(output_folder, "stat_table.tsv"), sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)

# Calls GIs and plots vs. null distribution
df$time1_fdr <- p.adjust(df$time1_pval, "BH", nrow(df))
df$time2_fdr <- p.adjust(df$time2_pval, "BH", nrow(df))
df$time1_is_gi <- df$time1_fdr < 0.1
df$time2_is_gi <- df$time2_fdr < 0.1
ggplot(df, aes(x = time1_expected_logFC, y = time1_observed_logFC)) +
  geom_point(aes(color = factor(time1_is_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Early timepoint expected log FC") +
  ylab("Early timepoint observed log FC") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_early_orientation_gis.png")))
ggplot(df, aes(x = time2_expected_logFC, y = time2_observed_logFC)) +
  geom_point(aes(color = factor(time2_is_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Late timepoint expected log FC") +
  ylab("Late timepoint observed log FC") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_late_orientation_gis.png")))

# Writes data to file
write.table(df, file.path(output_folder, paste0(which_cell_line, "_orientation_gi_calls.tsv")), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
saveRDS(time1_expected_sets, file.path(output_folder, paste0(which_cell_line, "_early_expected.rds")))
saveRDS(time1_observed_sets, file.path(output_folder, paste0(which_cell_line, "_early_observed.rds")))
saveRDS(time2_expected_sets, file.path(output_folder, paste0(which_cell_line, "_late_expected.rds")))
saveRDS(time2_observed_sets, file.path(output_folder, paste0(which_cell_line, "_late_observed.rds")))

# Gets gene-level GIs by checking which GIs agree
n <- nrow(df)/2
gene_df <- data.frame(gene1 = rep(NA, n), gene2 = rep(NA, n),
                      time1_expected_logFC1 = rep(NA, n), time1_expected_logFC2 = rep(NA, n),
                      time1_observed_logFC1 = rep(NA, n), time1_observed_logFC2 = rep(NA, n),
                      time1_fdr1 = rep(NA, n), time1_fdr2 = rep(NA, n),
                      time1_gi1 = rep(NA, n), time1_gi2 = rep(NA, n),
                      time2_expected_logFC1 = rep(NA, n), time2_expected_logFC2 = rep(NA, n),
                      time2_observed_logFC1 = rep(NA, n), time2_observed_logFC2 = rep(NA, n),
                      time2_fdr1 = rep(NA, n), time2_fdr2 = rep(NA, n),
                      time2_gi1 = rep(NA, n), time2_gi2 = rep(NA, n))
stat_table <- data.frame(gene1 = NA, gene2 = NA, 
                         n_expected_1 = NA, n_expected_2 = NA,
                         n_observed_1 = NA, n_observed_2 = NA)
counter <- 1
for (i in 1:nrow(df)) {
  gene1 <- df$gene1[i]
  gene2 <- df$gene2[i]
  temp <- gene_df[gene_df$gene1 == gene2 & gene_df$gene2 == gene1,]
  if (all(is.na(temp))) {
    row1 <- df[df$gene1 == gene1 & df$gene2 == gene2,]
    row2 <- df[df$gene2 == gene1 & df$gene1 == gene2,]
    expected1 <- as.numeric(row1$time1_expected_logFC)
    expected2 <- as.numeric(row2$time1_expected_logFC)
    observed1 <- as.numeric(row1$time1_observed_logFC)
    observed2 <- as.numeric(row2$time1_observed_logFC)
    fdr1 <- as.numeric(row1$time1_fdr)
    fdr2 <- as.numeric(row2$time1_fdr)
    gi1 <- row1$time1_is_gi
    gi2 <- row2$time1_is_gi
    time1_vec <- c(gene1, gene2, expected1, expected2, observed1, observed2, fdr1, fdr2, gi1, gi2)
    expected1 <- as.numeric(row1$time2_expected_logFC)
    expected2 <- as.numeric(row2$time2_expected_logFC)
    observed1 <- as.numeric(row1$time2_observed_logFC)
    observed2 <- as.numeric(row2$time2_observed_logFC)
    fdr1 <- as.numeric(row1$time2_fdr)
    fdr2 <- as.numeric(row2$time2_fdr)
    gi1 <- row1$time2_is_gi
    gi2 <- row2$time2_is_gi
    time2_vec <- c(expected1, expected2, observed1, observed2, fdr1, fdr2, gi1, gi2)
    gene_df[counter,] <- c(time1_vec, time2_vec)
    counter <- counter + 1
  }
}
gene_df[,3:8] <- sapply(gene_df[,3:8], as.numeric)
gene_df[,9:10] <- sapply(gene_df[,9:10], as.logical)
gene_df[,11:16] <- sapply(gene_df[,11:16], as.numeric)
gene_df[,17] <- as.logical(as.integer(gene_df[,17]))
gene_df[,18] <- as.logical(as.integer(gene_df[,18]))
gene_df$time1_mean_expected_logFC <- rowMeans(gene_df[,c("time1_expected_logFC1", "time1_expected_logFC2")])
gene_df$time1_mean_observed_logFC <- rowMeans(gene_df[,c("time1_observed_logFC1", "time1_observed_logFC2")])
gene_df$time2_mean_expected_logFC <- rowMeans(gene_df[,c("time2_expected_logFC1", "time2_expected_logFC2")])
gene_df$time2_mean_observed_logFC <- rowMeans(gene_df[,c("time2_observed_logFC1", "time2_observed_logFC2")])
gene_df$time1_both_gi <- gene_df$time1_gi1 & gene_df$time1_gi2
gene_df$time2_both_gi <- gene_df$time2_gi1 & gene_df$time2_gi2

# Flags GIs with opposite signs
time1_sign1 <- sign(gene_df$time1_observed_logFC1 - gene_df$time1_expected_logFC1)
time1_sign2 <- sign(gene_df$time1_observed_logFC2 - gene_df$time1_expected_logFC2)
time2_sign1 <- sign(gene_df$time2_observed_logFC1 - gene_df$time2_expected_logFC1)
time2_sign2 <- sign(gene_df$time2_observed_logFC2 - gene_df$time2_expected_logFC2)
time1_flagged <- time1_sign1 != time1_sign2 & gene_df$time1_both_gi
time2_flagged <- time2_sign1 != time2_sign2 & gene_df$time2_both_gi
gene_df$time1_both_gi[time1_flagged] <- FALSE
gene_df$time2_both_gi[time2_flagged] <- FALSE

# Plots gene-level GIs
ggplot(gene_df, aes(x = time1_mean_expected_logFC, y = time1_mean_observed_logFC)) +
  geom_point(aes(color = factor(time1_both_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Early expected log FC") +
  ylab("Early observed log FC") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_early_gene_gis.png")))
ggplot(gene_df, aes(x = time2_mean_expected_logFC, y = time2_mean_observed_logFC)) +
  geom_point(aes(color = factor(time2_both_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Late expected log FC") +
  ylab("Late observed log FC") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_late_gene_gis.png")))

# Writes data to file
write.table(gene_df, file.path(output_folder, paste0(which_cell_line, "_gene_gi_calls.tsv")), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Plots gene-level differences between orientation-level values for the early time
ggplot(gene_df, aes(x = time1_expected_logFC1, y = time1_expected_logFC2)) +
  geom_point(aes(color = factor(time1_both_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Early gene-level expected log FC (orientation 1)") +
  ylab("Early gene-level expected log FC (orientation 2)") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_early_orientation_expected_diffs.png")))
ggplot(gene_df, aes(x = time1_observed_logFC1, y = time1_observed_logFC2)) +
  geom_point(aes(color = factor(time1_both_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Early gene-level observed log FC (orientation 1)") +
  ylab("Early gene-level observed log FC (orientation 2)") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_early_orientation_observed_diffs.png")))
ggplot(gene_df, aes(x = time1_fdr1, y = time1_fdr2)) +
  geom_point(aes(color = factor(time1_both_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Early gene-level FDR (orientation 1)") +
  ylab("Early gene-level FDR (orientation 2)") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14)
ggsave(file.path(output_folder, paste0(which_cell_line, "_early_orientation_fdr_diffs.png")))
ggplot(gene_df, aes(x = time1_observed_logFC1 - time1_expected_logFC1, y = time1_observed_logFC2 - time1_expected_logFC2)) +
  geom_point(aes(color = factor(time1_both_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  geom_hline(yintercept = 0, linetype = 2, color = "Gray") +
  geom_vline(xintercept = 0, linetype = 2, color = "Gray") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Early gene-level residuals (orientation 1)") +
  ylab("Early gene-level residuals (orientation 2)") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_early_orientation_residual_diffs_both.png")))
ggplot(gene_df, aes(x = time1_observed_logFC1 - time1_expected_logFC1, y = time1_observed_logFC2 - time1_expected_logFC2)) +
  geom_point(aes(color = factor(time1_gi1 | time1_gi2, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  geom_hline(yintercept = 0, linetype = 2, color = "Gray") +
  geom_vline(xintercept = 0, linetype = 2, color = "Gray") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Early gene-level residuals (orientation 1)") +
  ylab("Early gene-level residuals (orientation 2)") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_early_orientation_residual_diffs_either.png")))

# Plots gene-level differences between orientation-level values for the late time
ggplot(gene_df, aes(x = time2_expected_logFC1, y = time2_expected_logFC2)) +
  geom_point(aes(color = factor(time2_both_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Late gene-level expected log FC (orientation 1)") +
  ylab("Late gene-level expected log FC (orientation 2)") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_late_orientation_expected_diffs.png")))
ggplot(gene_df, aes(x = time2_observed_logFC1, y = time2_observed_logFC2)) +
  geom_point(aes(color = factor(time2_both_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Late gene-level observed log FC (orientation 1)") +
  ylab("Late gene-level observed log FC (orientation 2)") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_late_orientation_observed_diffs.png")))
ggplot(gene_df, aes(x = time2_fdr1, y = time2_fdr2)) +
  geom_point(aes(color = factor(time2_both_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Late gene-level FDR (orientation 1)") +
  ylab("Late gene-level FDR (orientation 2)") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14)
ggsave(file.path(output_folder, paste0(which_cell_line, "_late_orientation_fdr_diffs.png")))
ggplot(gene_df, aes(x = time2_observed_logFC1 - time2_expected_logFC1, y = time2_observed_logFC2 - time2_expected_logFC2)) +
  geom_point(aes(color = factor(time2_both_gi, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  geom_hline(yintercept = 0, linetype = 2, color = "Gray") +
  geom_vline(xintercept = 0, linetype = 2, color = "Gray") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Late gene-level residuals (orientation 1)") +
  ylab("Late gene-level residuals (orientation 2)") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_late_orientation_residual_diffs_both.png")))
ggplot(gene_df, aes(x = time2_observed_logFC1 - time2_expected_logFC1, y = time2_observed_logFC2 - time2_expected_logFC2)) +
  geom_point(aes(color = factor(time2_gi1 | time2_gi2, labels = c("fdr > 0.1", "fdr < 0.1")))) +
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  geom_hline(yintercept = 0, linetype = 2, color = "Gray") +
  geom_vline(xintercept = 0, linetype = 2, color = "Gray") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Late gene-level residuals (orientation 1)") +
  ylab("Late gene-level residuals (orientation 2)") +
  labs(color = "Significant GI") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.9))
ggsave(file.path(output_folder, paste0(which_cell_line, "_late_orientation_residual_diffs_either.png")))

# Prints summary stats about orientation agreement
agreement_cor <- round(cor(gene_df$time1_expected_logFC1, gene_df$time1_expected_logFC2), 3)
cat(paste0("Correlation between each orientation's expected log FCs per gene (early): ", agreement_cor, "\n"))
agreement_cor <- round(cor(gene_df$time1_observed_logFC1, gene_df$time1_observed_logFC2), 3)
cat(paste0("Correlation between each orientation's observed log FCs per gene (early): ", agreement_cor, "\n"))
agreement_cor <- round(cor(gene_df$time2_expected_logFC1, gene_df$time2_expected_logFC2), 3)
cat(paste0("Correlation between each orientation's expected log FCs per gene (late): ", agreement_cor, "\n"))
agreement_cor <- round(cor(gene_df$time2_observed_logFC1, gene_df$time2_observed_logFC2), 3)
cat(paste0("Correlation between each orientation's observed log FCs per gene (late): ", agreement_cor, "\n"))
sink()

