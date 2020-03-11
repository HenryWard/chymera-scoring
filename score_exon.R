#!/usr/bin/env Rscript

library(ggplot2)
library(ggthemes)
library(dplyr)
library(scales)

######
# PARAMETER SETTING
######

# MODIFY THE TWO VARIABLES BELOW TO RUN THIS SCRIPT ON YOUR MACHINE!!!

# Sets working directory, input file location, output folder path, and cell line
setwd("/project/csbio/henry/Documents/projects/dual_guide/chymera_github")
output_folder <- file.path("output", "exon_deletion")

######
# MAIN SCRIPT
######

# Makes output folder if it doesn't exist
if(!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Sets important parameters
output_file <- file.path(output_folder, "results.txt")

# Sinks results to a text file
sink(output_file)

# Reads in and subsets data
input_file <- file.path("input", "exon_deletion_norm_counts.txt")
exon <- read.csv(input_file, sep = "\t", stringsAsFactors = FALSE)
exonic <- exon[exon$Cas9_Guide_Type == "exonic",]
exon <- exon[!(exon$Data_Subset %in% c("ExonDeletion_TKOv3", "ExonDeletion_Cas9Targeting")),]
intergenic <- exon[exon$Cas9_Guide_Type == "intergenic" & exon$Cpf1_Guide_Type == "intergenic",]
exon <- exon[!(exon$Cas9_Guide_Type == "intergenic" & exon$Cpf1_Guide_Type == "intergenic"),]
singles <- exon[exon$Data_Subset == "ExonDeletion_Ctrl",]
exon <- exon[exon$Data_Subset == "ExonDeletion",]

# Joins guides into one df
to_keep <- c("ID", "Event", "Target_gene", "Cas9_Guide_Type", "Cpf1_Guide_Type", "Cas9_Target_Site", "Cpf1_Target_Site", 
             "HAP1_logFC", "RPE1_logFC")
exon <- exon[,colnames(exon) %in% to_keep]
exon$type <- "exon"
exonic <- exonic[,colnames(exonic) %in% to_keep]
exonic$type <- "exonic"
singles <- singles[,colnames(singles) %in% to_keep]
singles$type <- "single"
intergenic <- intergenic[,colnames(intergenic) %in% to_keep]
intergenic$type <- "intergenic"
df <- rbind(exon, exonic, singles, intergenic)

# Compares two distributions
print(p.adjust(wilcox.test(exon$HAP1_logFC, intergenic$HAP1_logFC)$p.value), 3)
print(p.adjust(wilcox.test(singles$HAP1_logFC, intergenic$HAP1_logFC)$p.value), 3)
print(p.adjust(wilcox.test(exonic$HAP1_logFC, intergenic$HAP1_logFC)$p.value), 3)
print(p.adjust(wilcox.test(exon$RPE1_logFC, intergenic$RPE1_logFC)$p.value), 3)
print(p.adjust(wilcox.test(singles$RPE1_logFC, intergenic$RPE1_logFC)$p.value), 3)
print(p.adjust(wilcox.test(exonic$RPE1_logFC, intergenic$RPE1_logFC)$p.value), 3)

# Gets target-site-specific intronic-intergenic scores 
events <- 
  singles %>% 
  group_by(Event, Cas9_Target_Site, Cpf1_Target_Site) %>%
  summarise(mean_HAP1 = mean(HAP1_logFC),
            mean_RPE1 = mean(RPE1_logFC),
            sd_HAP1 = sd(HAP1_logFC),
            sd_RPE1 = sd(RPE1_logFC),
            total = n())

# Gets frequency of observations for intergenic log FCs
hap1_freq <- data.frame(val = round(seq(-7.3, 5.8, 0.01), 2),
                        num = rep(0, length(seq(-7.3, 5.8, 0.01))))
rpe1_freq <- data.frame(val = round(seq(-7.3, 5.8, 0.01), 2),
                        num = rep(0, length(seq(-7.3, 5.8, 0.01))))
hap1_null <- round(intergenic$HAP1_logFC[order(intergenic$HAP1_logFC)], 2)
rpe1_null <- round(intergenic$RPE1_logFC[order(intergenic$RPE1_logFC)], 2)
for (i in 1:length(hap1_null)) {
  ind <- hap1_freq$val %in% hap1_null[i]
  hap1_freq$num[ind] <- hap1_freq$num[ind] + 1
}
for (i in 1:length(rpe1_null)) {
  ind <- rpe1_freq$val %in% rpe1_null[i]
  rpe1_freq$num[ind] <- rpe1_freq$num[ind] + 1
}
hap1_freq$neg_FDR <- cumsum(hap1_freq$num) / sum(hap1_freq$num)
hap1_freq$pos_FDR <- rev(cumsum(rev(hap1_freq$num)) / sum(hap1_freq$num))
rpe1_freq$neg_FDR <- cumsum(rpe1_freq$num) / sum(rpe1_freq$num)
rpe1_freq$pos_FDR <- rev(cumsum(rev(rpe1_freq$num)) / sum(rpe1_freq$num))

# Gets frequency of observations for intronic-intronic guide data relative to intergenic controls
# and relative to intronic-intergenic guides
exon$hap1_fdr <- NA
exon$rpe1_fdr <- NA
exon$hap1_cas9_single <- NA
exon$hap1_cpf1_single <- NA
exon$rpe1_cas9_single <- NA
exon$rpe1_cpf1_single <- NA
exon$cas9_n_single <- NA
exon$cpf1_n_single <- NA
for (i in 1:nrow(exon)) {
  hap1_subset <- hap1_freq[hap1_freq$val == round(exon$HAP1_logFC[i], 2),]
  rpe1_subset <- rpe1_freq[rpe1_freq$val == round(exon$RPE1_logFC[i], 2),]
  exon$hap1_fdr[i] <- min(hap1_subset$neg_FDR, hap1_subset$pos_FDR) 
  exon$rpe1_fdr[i] <- min(rpe1_subset$neg_FDR, rpe1_subset$pos_FDR)
  cas9_set <- events[events$Cas9_Target_Site == exon$Cas9_Target_Site[i],]
  cpf1_set <- events[events$Cpf1_Target_Site == exon$Cpf1_Target_Site[i],]
  exon$hap1_cas9_single[i] <- mean(cas9_set$mean_HAP1)
  exon$hap1_cpf1_single[i] <- mean(cpf1_set$mean_HAP1)
  exon$rpe1_cas9_single[i] <- mean(cas9_set$mean_RPE1)
  exon$rpe1_cpf1_single[i] <- mean(cpf1_set$mean_RPE1)
  exon$cas9_n_single <- nrow(cas9_set)
  exon$cpf1_n_single <- nrow(cpf1_set)
}
exon$hap1_sig <- exon$hap1_fdr < 0.025
exon$rpe1_sig <- exon$rpe1_fdr < 0.025
exon$hap1_mean_single <- rowMeans(exon[,(colnames(exon) %in% c("hap1_cas9_single", "hap1_cpf1_single"))])
exon$rpe1_mean_single <- rowMeans(exon[,(colnames(exon) %in% c("rpe1_cas9_single", "rpe1_cpf1_single"))])
exon$hap1_residual_single <- exon$HAP1_logFC - exon$hap1_mean_single
exon$rpe1_residual_single <- exon$RPE1_logFC - exon$rpe1_mean_single

# Gets frequency of observations for exonic-intergenic guide data relative to intergenic controls
exonic$hap1_fdr <- NA
exonic$rpe1_fdr <- NA
for (i in 1:nrow(exonic)) {
  hap1_subset <- hap1_freq[hap1_freq$val == round(exonic$HAP1_logFC[i], 2),]
  rpe1_subset <- rpe1_freq[rpe1_freq$val == round(exonic$RPE1_logFC[i], 2),]
  exonic$hap1_fdr[i] <- min(hap1_subset$neg_FDR, hap1_subset$pos_FDR) 
  exonic$rpe1_fdr[i] <- min(rpe1_subset$neg_FDR, rpe1_subset$pos_FDR) 
}
exonic$hap1_sig <- exonic$hap1_fdr < 0.025
exonic$rpe1_sig <- exonic$rpe1_fdr < 0.025

# Gets frequency of observations for intronic-intergenic guide data relative to intergenic controls
singles$hap1_fdr <- NA
singles$rpe1_fdr <- NA
for (i in 1:nrow(singles)) {
  hap1_subset <- hap1_freq[hap1_freq$val == round(singles$HAP1_logFC[i], 2),]
  rpe1_subset <- rpe1_freq[rpe1_freq$val == round(singles$RPE1_logFC[i], 2),]
  singles$hap1_fdr[i] <- min(hap1_subset$neg_FDR, hap1_subset$pos_FDR) 
  singles$rpe1_fdr[i] <- min(rpe1_subset$neg_FDR, rpe1_subset$pos_FDR) 
}
singles$hap1_sig <- singles$hap1_fdr < 0.025
singles$rpe1_sig <- singles$rpe1_fdr < 0.025

# Gets boundaries of significant guides
exon_hap1_min <- max(exon$HAP1_logFC[exon$hap1_sig & exon$HAP1_logFC < 0])
exon_hap1_max <- min(exon$HAP1_logFC[exon$hap1_sig & exon$HAP1_logFC > 0])
exon_rpe1_min <- max(exon$RPE1_logFC[exon$rpe1_sig & exon$RPE1_logFC < 0])
exon_rpe1_max <- min(exon$RPE1_logFC[exon$rpe1_sig & exon$RPE1_logFC > 0])
exonic_hap1_min <- max(exonic$HAP1_logFC[exonic$hap1_sig & exonic$HAP1_logFC < 0])
exonic_hap1_max <- min(exonic$HAP1_logFC[exonic$hap1_sig & exonic$HAP1_logFC > 0])
exonic_rpe1_min <- max(exonic$RPE1_logFC[exonic$rpe1_sig & exonic$RPE1_logFC < 0])
exonic_rpe1_max <- min(exonic$RPE1_logFC[exonic$rpe1_sig & exonic$RPE1_logFC > 0])
single_hap1_min <- max(singles$HAP1_logFC[singles$hap1_sig & singles$HAP1_logFC < 0])
single_hap1_max <- min(singles$HAP1_logFC[singles$hap1_sig & singles$HAP1_logFC > 0])
single_rpe1_min <- max(singles$RPE1_logFC[singles$rpe1_sig & singles$RPE1_logFC < 0])
single_rpe1_max <- min(singles$RPE1_logFC[singles$rpe1_sig & singles$RPE1_logFC > 0])

print(exon_hap1_min)
print(exon_hap1_max)
print(exon_rpe1_min)
print(exon_rpe1_max)

# Plots distributions of intronic-intronic vs. intronic-intergenic vs. intergenic-intergenic data
pal <- hue_pal()(4)
ggplot(df, aes(HAP1_logFC)) +
  geom_density(aes(color = factor(type, labels = 
                                    c("Intronic-intronic", "Exonic-intergenic", "Intergenic-intergenic",  "Intronic-intergenic")))) +
  geom_vline(xintercept = exon_hap1_min, color = pal[1], linetype = 2) +
  geom_vline(xintercept = exon_hap1_max, color = pal[1], linetype = 2) +
  geom_vline(xintercept = exonic_hap1_min, color = pal[2], linetype = 2) +
  geom_vline(xintercept = exonic_hap1_max, color = pal[2], linetype = 2) +
  geom_vline(xintercept = single_hap1_min, color = pal[3], linetype = 2) +
  geom_vline(xintercept = single_hap1_max, color = pal[3], linetype = 2) +
  xlab("HAP1 logFC") +
  ylab("Density") +
  labs(color = "Guide type") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.85, 0.85))
ggsave(file.path(output_folder, "hap1_distributions.png"))
ggplot(df, aes(RPE1_logFC)) +
  geom_density(aes(color = factor(type, labels = 
                                    c("Intronic-intronic", "Exonic-intergenic", "Intergenic-intergenic", "Intronic-intergenic")))) +
  geom_vline(xintercept = exon_rpe1_min, color = pal[1], linetype = 2) +
  geom_vline(xintercept = exon_rpe1_max, color = pal[1], linetype = 2) +
  geom_vline(xintercept = exonic_rpe1_min, color = pal[2], linetype = 2) +
  geom_vline(xintercept = exonic_rpe1_max, color = pal[2], linetype = 2) +
  geom_vline(xintercept = single_rpe1_min, color = pal[3], linetype = 2) +
  geom_vline(xintercept = single_rpe1_max, color = pal[3], linetype = 2) +
  xlab("RPE1 logFC") +
  ylab("Density") +
  labs(color = "Guide type") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.85, 0.85))
ggsave(file.path(output_folder, "rpe1_distributions.png"))

# Plots intronic-intronic vs. intronic-intergenic log FCs
exon$orientation <- "Cas9-up, Cpf1-down"
exon$orientation[exon$Cas9_Guide_Type == "intronic_DN"] <- "Cas9-down, Cpf1-up"
ggplot(exon, aes(x = hap1_mean_single, y = HAP1_logFC)) +
  geom_point(aes(color =  factor(hap1_sig, labels = c("fdr > 0.05", "fdr < 0.05"))), alpha = 0.8) +
  geom_smooth(method = "lm") +
  #geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Intronic-intergenic HAP1 log FC") +
  ylab("Intronic-intronic HAP1 log FC") +
  labs(color = "Significant Guide Pair") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.85))
ggsave(file.path(output_folder, "hap1_single_comparison.png"))
ggplot(exon, aes(x = rpe1_mean_single, y = RPE1_logFC)) +
  geom_point(aes(color =  factor(rpe1_sig, labels = c("fdr > 0.05", "fdr < 0.05"))), alpha = 0.8) +
  geom_smooth(method = "lm") +
  #geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Intronic-intergenic RPE1 log FC") +
  ylab("Intronic-intronic RPE1 log FC") +
  labs(color = "Significant Guide Pair") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.85))
ggsave(file.path(output_folder, "rpe1_single_comparison.png"))

# Gets some stats
print(sd(exon$hap1_mean_single))
print(sd(exon$HAP1_logFC))
print(sd(exon$rpe1_mean_single))
print(sd(exon$RPE1_logFC))
print(cor(exon$hap1_mean_single, exon$HAP1_logFC))
print(cor(exon$rpe1_mean_single, exon$RPE1_logFC))
print(cor(exon$hap1_mean_single, exon$rpe1_mean_single))
print(cor(exon$HAP1_logFC, exon$RPE1_logFC))

exon$exonic_HAP1_logFC <- NA
exon$exonic_RPE1_logFC <- NA
for (i in 1:nrow(exon)) {
  exon$exonic_HAP1_logFC[i] <- mean(exonic$HAP1_logFC[exonic$Event == exon$Event[i]], na.rm = TRUE)
  exon$exonic_RPE1_logFC[i] <- mean(exonic$RPE1_logFC[exonic$Event == exon$Event[i]], na.rm = TRUE)
}

ggplot(exon, aes(x = exonic_HAP1_logFC, y = HAP1_logFC)) +
  geom_point(aes(color =  factor(hap1_sig, labels = c("fdr > 0.05", "fdr < 0.05"))), alpha = 0.8) +
  geom_smooth(method = "lm") +
  #geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Exonic-intergenic HAP1 log FC") +
  ylab("Intronic-intronic HAP1 log FC") +
  labs(color = "Significant Guide Pair") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.85))
ggsave(file.path(output_folder, "hap1_exonic_comparison.png"))
ggplot(exon, aes(x = exonic_RPE1_logFC, y = RPE1_logFC)) +
  geom_point(aes(color =  factor(rpe1_sig, labels = c("fdr > 0.05", "fdr < 0.05"))), alpha = 0.8) +
  geom_smooth(method = "lm") +
  #geom_abline(slope = 1, intercept = 0, size = 1.5, color = "Blue") +
  scale_color_manual(values = c("Gray", "Red")) +
  xlab("Exonic-intergenic RPE1 log FC") +
  ylab("Intronic-intronic RPE1 log FC") +
  labs(color = "Significant Guide Pair") +
  theme_tufte(base_size = 14) +
  theme(legend.position = c(0.15, 0.85))
ggsave(file.path(output_folder, "rpe1_exonic_comparison.png"))

# Writes data to file and ends sinking
write.table(rbind(exonic, singles), file.path(output_folder, "other_guides.tsv"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(exon, file.path(output_folder, "intronic_intronic_guides.tsv"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Ends sinking
sink()