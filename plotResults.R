library(ggplot2)


data <- read.csv("data/resultsBenchmark.csv", sep = ';')

colnames(data)[1] <- "lfc"
ggplot(data, aes(color = method, x = factor(fdr), y = precision)) + 
  geom_point(size = 2) +
  facet_wrap(~lfc + density)


ggplot(data[data$time > 0,], aes(x = genes, color = factor(density), y = time)) + geom_point(size = 5)
