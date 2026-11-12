library(dplyr)
library(ggplot2)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(readr)
library(gtools)
library(cowplot)

#working directory: mutation-spectrum-misspecification-in-dNdS-estimates-and-the-extended-Nei-Gojobori-method/
#This script is responsible for:
# - converting summary_of_mutation_spectrum_and_omega_ratio.csv to summary_of_mutation_spectrum_and_omega_ratio_updated.csv
# - Figure making for Figure 2, 3, 4, S1, S2
# - Table S2

#data loading
data = read.csv("data/summary_of_mutation_spectrum_and_omega_ratio.csv")

#discard annotation rows and some species
data = data[-(47:48),]
data = data[-37,]
data = data[-27,]
data = data[-17,]
data = data[-11,]


###############################
#data processing
colnames(data) <- c("Groups", "Species", "ATtoTA", "ATtoCG","ATtoGC",
                    "CGtoGC", "CGtoAT", "CGtoTA", "GC_AT_mutation_rate_ratio", "omega_ratio")

data$inverse_of_omega_ratio = 1/data$omega_ratio

data$sd1 = apply(data[,3:5], 1, sd)
data$sd2 = apply(data[,6:8], 1, sd)
data$mutation_direction_bias = 3*(0.5*data$sd1 + 0.5*data$sd2)
summary(data$mutation_direction_bias)

data$omega_fold_difference[data$omega_ratio>1] <- data$omega_ratio[data$omega_ratio>1]
data$omega_fold_difference[data$omega_ratio<1] <- data$inverse_of_omega_ratio[data$omega_ratio<1]
summary(data$omega_fold_difference)

rownames(data) <- 1:42
data$Groups = factor(data$Groups, levels = c("Bacteria", "Unicellular eukaryotes", "Multicellular eukaryotes"))

data$estimate[data$omega_ratio>1] <- "Overestimate"
data$estimate[data$omega_ratio<1] <- "Underestimate"

bacteria <- data[data$Groups == "Bacteria",]
unicellular <- data[data$Groups == "Unicellular eukaryotes",]
multicellular <- data[data$Groups == "Multicellular eukaryotes",]

table(data$estimate)
table(bacteria$estimate)
table(unicellular$estimate)
table(multicellular$estimate)


data$species_name = c("Agrobacterium tumefaciens",
                      "Arthrobacter sp.",
                      "Bacillus subtilis",
                      "Burkhoderia cenocepacia",
                      "Caulobacter crescentus",
                      "Colwellia psychrerythraea",
                      "Deinococcus radiodurans",
                      "Escherichia coli",
                      "Flavobacterium sp.",
                      "Gemmata obscuriglobus",
                      "Kineococcus radiotolerans",
                      "Lactococcus lactis",
                      "Mesoplasma florum",
                      "Micrococcus sp.",
                      "Mycobacterium smegmatis",
                      "Ruegeria pomeroyi",
                      "Staphylococcus aureus",
                      "Staphylococcus epidermidis",
                      "Vibrio cholerae",
                      "Vibrio fischeri",
                      "JCVI-syn3B",
                      "JCVI-syn1.0",
                      "Arabidopsis thaliana",
                      "Caenorhabditis elegans",
                      "Drosophila melanogaster",
                      "Ficedula albicollis",
                      "Homo sapiens",
                      "Mus musculus",
                      "Oryza sativa",
                      "Pan troglodytes",
                      "Chlamydomonas reinhardtii",
                      "Plasmodium falciparum",
                      "Rhodotorula toruloides",
                      "Schizosaccharomyces pombe",
                      "Sphaeroforma arctica",
                      "Saccharomyces cerevisiae (Cu)",
                      "Saccharomyces cerevisiae (Li)",
                      "Saccharomyces cerevisiae (Na)",
                      "Saccharomyces cerevisiae (YNB)",
                      "Saccharomyces cerevisiae (YPD)",
                      "Saccharomyces cerevisiae (YPL)",
                      "Saccharomyces cerevisiae (YPX)"
)

###############################
#Figure making for Figure 3: omega ratio and omega fold change

purple = "#9896C2"
yellow = "#EBBD61"
green = "#7FB361"
colors = c(purple, orange, green2)

#Figure3A: omega ratio
p <- ggplot(data, aes(x = omega_ratio, fill = Groups)) +
  geom_histogram(binwidth = 0.05, boundary=0.5, alpha = 0.9)+
  coord_trans(x = "log")+
  theme_minimal()+
  scale_fill_manual(values = colors)+
  scale_x_continuous(breaks = seq(0.5, 3, by = 0.5),
                     limits = c(0.5, 3))+
  scale_y_continuous(breaks = seq(0, 10, by = 2),
                     limits = c(0, 9))+
  labs(x = expression(italic("\U03C9")~"ratio"), 
       y = 'Count')+
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.x  = element_text(size = 11, color = "black"),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, color = "black"),
  )

grDevices::quartz(
  file = "figures/omega_ratio.pdf",
  type = "pdf",
  width = 6,
  height = 3
)
print(p)
dev.off()

#Figure 3B: omega fold diff
p <- ggplot(data, aes(x = omega_fold_difference, fill = Groups)) +
  geom_histogram(binwidth = 0.1, boundary=1, alpha = 0.9)+
  theme_minimal()+
  scale_fill_manual(values = colors)+
  scale_x_continuous()+
  scale_y_continuous(breaks = seq(0, 12, by = 2),
                     limits = c(0, 12))+
  labs(x = expression(italic("\U03C9")~"fold difference"), 
       y = 'Count')+
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.x  = element_text(size = 11, color = "black"),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, color = "black"),
  )

grDevices::quartz(
  file = "figures/omega_fold_difference.pdf",
  type = "pdf",
  width = 6,
  height = 3
)
print(p)
dev.off()

###############################
#Figure making for Figure 2: mutation direction bias and GC/AT mutation rate ratio

#Fig 2A
p <- ggplot(data, aes(x = mutation_direction_bias, fill = Groups)) +
  geom_histogram(binwidth = 0.05, boundary = c(0.35), alpha = 0.9)+
  theme_minimal()+
  scale_fill_manual(values = colors1)+
  scale_x_continuous(breaks = seq(0.2, 1.4, by = 0.2))+
  labs(x = ' Mutation direction bias', y = 'Count')+
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.x  = element_text(size = 13, color = "black"),
    axis.text.y  = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 13, color = "black")
  )

grDevices::quartz(
  file = "figures/mutation_direction_bias.pdf",
  type = "pdf",
  width = 7,
  height = 4
)
print(p)
dev.off()


#Fig 2B
p <- ggplot(data, aes(x = GC_AT_mutation_rate_ratio, fill = Groups)) +
  geom_histogram(binwidth = 1, boundary = 0, alpha = 0.9)+
  theme_minimal()+
  scale_fill_manual(values = colors1)+
  scale_y_continuous(breaks = seq(0, 14, by = 2))+
  scale_x_continuous(breaks = seq(0, 45, by = 5),
                     limits = c(0, 45))+
  labs(x = 'GC/AT mutation rate ratio', y = 'Count')+
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.x  = element_text(size = 13, color = "black"),
    axis.text.y  = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    legend.position = "none"
  )

grDevices::quartz(
  file = "figures/GC_AT_mutation_rate_ratio.pdf",
  type = "pdf",
  width = 7,
  height = 4
)
print(p)
dev.off()

###############################
#Table making for Table S2: Statistics
statistics_omega_ratio <- data %>% group_by(Groups) %>% summarise(
  min = min(omega_ratio),
  one_fourth = quantile(omega_ratio, 0.25),
  mean = mean(omega_ratio),
  median = median(omega_ratio),
  three_fourth = quantile(omega_ratio, 0.75),
  max = max(omega_ratio),
  sd = sd(omega_ratio),
)


statistics_omega_fold_difference <- data %>% group_by(Groups) %>% summarise(
  min = min(omega_fold_difference),
  one_fourth = quantile(omega_fold_difference, 0.25),
  mean = mean(omega_fold_difference),
  median = median(omega_fold_difference),
  three_fourth = quantile(omega_fold_difference, 0.75),
  max = max(omega_fold_difference),
  sd = sd(omega_fold_difference),
)


statistics_total_omega_ratio <- data %>%
  summarise(
    Groups = "Total",
    min = min(omega_ratio),
    one_fourth = quantile(omega_ratio, 0.25),
    mean = mean(omega_ratio),
    median = median(omega_ratio),
    three_fourth = quantile(omega_ratio, 0.75),
    max = max(omega_ratio),
    sd = sd(omega_ratio)
  )

statistics_total_omega_fold_difference <- data %>%
  summarise(
    Groups = "Total",
    min = min(omega_fold_difference),
    one_fourth = quantile(omega_fold_difference, 0.25),
    mean = mean(omega_fold_difference),
    median = median(omega_fold_difference),
    three_fourth = quantile(omega_fold_difference, 0.75),
    max = max(omega_fold_difference),
    sd = sd(omega_fold_difference)
  )

statistics <- bind_rows(statistics_omega_ratio, statistics_omega_fold_difference, 
                        statistics_total_omega_ratio, statistics_total_omega_fold_difference)

write.csv(statistics, "data/statistics.csv")


###############################
#Figure making for Figure 4, Figure S1, Figure S2: omega ratio for individual genes

#data loading
folder = "data/omega_ratio_for_individual_genes"
files <- list.files(path = folder, full.names = TRUE)
files <- mixedsort(files)
data_list <- lapply(files, function(f){
  as.numeric(readLines(f))
})
max_len <- max(sapply(data_list, length))
data_list_padded <- lapply(data_list, function(x) {
  length(x) <- max_len
  x
})

summary_data <- as.data.frame(do.call(rbind, data_list_padded), stringsAsFactors = FALSE)
summary_data <- as.data.frame(t(summary_data))
summary_data <- as.data.frame(apply(summary_data, 2, as.numeric))

species_list = c("Agrobacterium tumefaciens",
                 "Arthrobacter",
                 "Bacillus subtilis",
                 "Burkhoderia cenocepacia",
                 "Caulobacter crescentus",
                 "Colwellia psychrerythraea",
                 "Deinococcus radiodurans",
                 "Escherichia coli",
                 "Flavobacterium",
                 "Gemmata obscuriglobus",
                 "Kineococcus radiotolerans",
                 "Lactococcus lactis",
                 "Mesoplasma florum",
                 "Micrococcus",
                 "Mycobacterium smegmatis",
                 "Ruegeria pomeroyi",
                 "Staphylococcus aureus",
                 "Staphylococcus epidermidis",
                 "Vibrio cholerae",
                 "Vibrio fischeri",
                 "",
                 "",
                 "Arabidopsis thaliana",
                 "Caenorhabditis elegans",
                 "Drosophila melanogaster",
                 "Ficedula albicollis",
                 "Homo sapiens",
                 "Mus musculus",
                 "Oryza sativa",
                 "Pan troglodytes",
                 "Chlamydomonas reinhardtii",
                 "Plasmodium falciparum",
                 "Rhodotorula toruloides",
                 "Schizosaccharomyces pombe",
                 "Sphaeroforma arctica",
                 "Saccharomyces cerevisiae",
                 "Saccharomyces cerevisiae",
                 "Saccharomyces cerevisiae",
                 "Saccharomyces cerevisiae",
                 "Saccharomyces cerevisiae",
                 "Saccharomyces cerevisiae",
                 "Saccharomyces cerevisiae"
                 )

environment <- c("", "sp.", rep("", 6), "sp.", rep("", 4), "sp.", rep("", 6),
                 "JCVI-syn3B", "JCVI-syn1.0", rep("", 13),
                "in CuSO4", "in LiCl", "in NaCl", "in YNB", "in YPD", "in YPL", "in YPX")


#Figure S2
summary_individual_genes <- as.data.frame(matrix(nrow=42, ncol=10))
colnames(summary_individual_genes) <- c("min", "Q1", "median", "mean", "Q3", "max", 
                                        "sd", "logsd", "trans_sd", "species")

for (i in 1:42){
  summary_individual_genes[i,1] <- min(summary_data[[i]], na.rm=TRUE)
  summary_individual_genes[i,2] <- quantile(summary_data[[i]], 0.25, na.rm=TRUE)
  summary_individual_genes[i,3] <- median(summary_data[[i]], na.rm=TRUE)
  summary_individual_genes[i,4] <- mean(summary_data[[i]], na.rm=TRUE)
  summary_individual_genes[i,5] <- quantile(summary_data[[i]], 0.75, na.rm=TRUE)
  summary_individual_genes[i,6] <- max(summary_data[[i]], na.rm=TRUE)
  summary_individual_genes[i,7] <- sd(summary_data[[i]], na.rm=TRUE)
  summary_individual_genes[i,8] <- sd(log(summary_data[[i]]), na.rm=TRUE)
  summary_individual_genes[i,9] <- sd(ifelse(summary_data[[i]] > 1,
                                             summary_data[[i]],
                           -1 / summary_data[[i]]), na.rm=TRUE)
}

summary_individual_genes$species <- data$species_name


data$sd_of_signed_omega_fold_difference <- summary_individual_genes$trans_sd

p <- ggplot(data, aes(x = sd_of_signed_omega_fold_difference, fill = Groups)) +
  geom_histogram(binwidth = 0.05, boundary=0, alpha = 0.9)+
  theme_minimal()+
  scale_fill_manual(values = colors1)+
  scale_x_continuous()+
  labs(x = expression("Standard deviation of the signed"~italic("\U03C9")~"fold difference across genes"), 
       y = 'Count')+
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.x  = element_text(size = 11, color = "black"),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, color = "black")
  )

grDevices::quartz(
  file = "figures/sd_of_signed_omega_fold_difference_across_genes.pdf",
  type = "pdf",
  width = 8,
  height = 5
)
print(p)
dev.off()


#Figure S1
plots <- list()

for (i in 1:42){
  g <- ggplot(summary_data, aes(x = .data[[colnames(summary_data)[i]]])) +
    geom_density(size = 1, color = 'black') + # 总体分布
    geom_vline(xintercept = data$omega_ratio[i], color = rgb(0.444,0.677,0.470), size = 0.7, alpha = 0.9, linetype = "longdash")+
    geom_vline(xintercept = 1, color = rgb(0.509,0.420,0.680), size = 0.7, alpha = 0.9, linetype = "longdash")+
    coord_trans(x = "log", xlim = c(0.09,4.02), ylim = c(0,20))+
    theme_minimal_grid() +
    scale_x_continuous(breaks = c(0.1,0.25,0.5,1,2,4))+
    labs(
      x = expression(italic("\U03C9")~"ratio"),
      y = "Probability density",
      title = bquote(italic(.(species_list[i]))~.(environment[i]))
         ) +
    theme(
      text = element_text(family = "Times New Roman"),
      axis.text.x  = element_text(size = 13, color = "black"),
      axis.text.y  = element_text(size = 13, color = "black"),
      axis.title.x = element_text(size = 17.5, color = "black"),
      axis.title.y = element_text(size = 17.5, color = "black"),
      plot.title = element_text(size = 17.5, color = "black", hjust = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    )
  plots[[i]] <- g
}

p <- plot_grid(plotlist = plots[1:15], ncol = 3)
grDevices::quartz(
  file = "figures/omega_ratio_for_individual_genes_1.pdf",
  type = "pdf",
  width = 13,
  height = 17.2
)
print(p)
dev.off()

p <- plot_grid(plotlist = plots[16:30], ncol = 3)
grDevices::quartz(
  file = "figures/omega_ratio_for_individual_genes_2.pdf",
  type = "pdf",
  width = 13,
  height = 17.2
)
print(p)
dev.off()

p <- plot_grid(plotlist = plots[31:42], ncol = 3)
grDevices::quartz(
  file = "figures/omega_ratio_for_individual_genes_3.pdf",
  type = "pdf",
  width = 13,
  height = 13.8
)
print(p)
dev.off()

# Figure 4: selected cases
selected_plots <- list()
k <- 0

for (i in c(39,8,12,23,27,13)){
  g <- ggplot(summary_data, aes(x = .data[[colnames(summary_data)[i]]])) +
    geom_density(size = 1, color = "black") + # 总体分布
    geom_vline(xintercept = data$omega_ratio[i], color = rgb(0.444,0.677,0.470), size = 0.7, alpha = 0.9, linetype = "longdash")+
    geom_vline(xintercept = 1, color = rgb(0.509,0.420,0.680), size = 0.7, alpha = 0.9, linetype = "longdash")+
    coord_trans(x = "log", xlim = c(0.2,4.02), ylim = c(0, 20))+
    theme_minimal_grid() +
    scale_x_continuous(breaks = c(0.1,0.25,0.5,1,2,4))+
    labs(
      x = expression(italic("\U03C9")~"ratio"),
      y = "Probability density",
      title = bquote(italic(.(species_list[i]))~.(environment[i]))
    )+
    theme(
      text = element_text(family = "Times New Roman"),
      axis.text.x  = element_text(size = 13, color = "black"),
      axis.text.y  = element_text(size = 13, color = "black"),
      axis.title.x = element_text(size = 17.5, color = "black"),
      axis.title.y = element_text(size = 17.5, color = "black"),
      plot.title = element_text(size = 17.5, color = "black", hjust = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    )
  k <- k + 1
  selected_plots[[k]] <- g
}

p <- plot_grid(plotlist = selected_plots, ncol = 2)

grDevices::quartz(
  file = "figures/omega_ratio_for_individual_genes_selected.pdf",
  type = "pdf",
  width = 8,
  height = 10
)
print(p)
dev.off()

write.csv(data, "data/summary_of_mutation_spectrum_and_omega_ratio_updated.csv", row.names = FALSE)
