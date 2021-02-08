library(patchwork)
library(ggplot2)
library(ggpubr)
library(stringr)


######### with thrid strategy ; 
# reading results and subsetting settings of interest


data$fdr <- paste("FDR :", data$fdr)
data$density <- paste("Density :", data$density)
data$Strategy <- str_replace(data$Strategy, "before_testing", "No tests")
data$Strategy <- str_replace(data$Strategy, "testing", "Testing")


connecTF <- ggplot(data, aes(color = Strategy, fill = Strategy, x = Strategy, y = precision, label = N_edges)) + 
  geom_boxplot(size = 0.55, alpha = 0.65, width = 0.25) + geom_jitter(size = 0.5, width = 0.05)+ 
  geom_text(check_overlap = TRUE, size =2.7, nudge_x = 0.37) +
  facet_wrap(~ density + factor(fdr), nrow = 2) + scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") + stat_compare_means(
    aes(x = Strategy, y = precision),
    comparisons = list(c("Testing", "No tests"), c("Testing", "same_edges")), method = "wilcox.test", paired = FALSE) + 
  ggtitle("A. thaliana : precision benchmarked on connecTF")  + 
  labs(caption = "Wilcox.test, N = 20") + xlab("") + ylab("Precision")+ theme(title = element_text(size = 13))

connecTF 

################## Two strategies as in the paper

# reading results and subsetting settings of interest
data <- read.csv("benchmark_ecoli_1000Trees_10N.csv")
data <- data[data$fdr <0.1,]
data <- data[data$density < 0.0075,]
data <- data[data$Strategy != "same n edges",]
data$fdr <- paste("FDR :", data$fdr)
data$density <- paste("Density :", data$density)
data$Strategy <- str_replace(data$Strategy, "before testing", "No tests")
data$Strategy <- str_replace(data$Strategy, "testing", "Testing")


ecoli <- ggplot(data, aes(color = Strategy, fill = Strategy, x = Strategy, y = precision, label = N_edges)) + 
  geom_boxplot(size = 0.5, alpha = 0.65, width = 0.25) + geom_jitter(size = 0.5, width = 0.05) +
  geom_text(check_overlap = TRUE, size =2.7, nudge_x = 0.37) +
  facet_wrap(~density + factor(fdr), nrow = 2) + scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") + stat_compare_means(
    aes(x = Strategy, y = precision),
    comparisons = list(c("Testing", "No tests")), method = "wilcox.test", paired = FALSE) + 
  ggtitle("E. coli : precision benchmarked on RegulonDB")  + theme(legend.position = "none", title = element_text(size = 13))+ 
  labs(caption = "Wilcox.test, N = 20") + xlab("") + ylab("Precision") 


data <- read.csv("benchmark_athaliana_1000Trees_CvsH_noDap.csv")
data <- data[data$Strategy != "same_edges",]
data <- data[data$lfc <2,]

# testing Strategy compared to the same networks before testing for all FDRs 
# (already delt with in eColi script but not A thaliana, that inferred new GENIE3s for different FDRs,
# but we want to compare testing to the same prior networks)

for(d in c(0.01,0.03)){
  no_tests <- data[data$density == d & data$fdr == 0.005 & data$Strategy == "before_testing",]$precision
  data[data$density == d & data$fdr == 0.01& data$Strategy == "before_testing", "precision"] <- no_tests
  data[data$density == d & data$fdr == 0.05& data$Strategy == "before_testing", "precision"] <- no_tests
}


data$fdr <- paste("FDR :", data$fdr)
data$density <- paste("Density :", data$density)
data$Strategy <- str_replace(data$Strategy, "before_testing", "No tests")
data$Strategy <- str_replace(data$Strategy, "testing", "Testing")


connecTF <- ggplot(data, aes(color = Strategy, fill = Strategy, x = Strategy, y = precision, label = N_edges)) + 
  geom_boxplot(size = 0.55, alpha = 0.65, width = 0.25) + geom_jitter(size = 0.5, width = 0.05)+ 
  geom_text(check_overlap = TRUE, size =2.7, nudge_x = 0.37) +
  facet_wrap(~ density + factor(fdr), nrow = 2) + scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") + stat_compare_means(
    aes(x = Strategy, y = precision),
    comparisons = list(c("Testing", "No tests")), method = "wilcox.test", paired = FALSE) + 
  ggtitle("A. thaliana : precision benchmarked on connecTF")  + 
  labs(caption = "Wilcox.test, N = 20") + xlab("") + ylab("Precision")+ theme(title = element_text(size = 13))

ecoli + connecTF 







