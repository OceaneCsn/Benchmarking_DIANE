library(DIANE)
data("abiotic_stresses")
data("gene_annotations")
data("regulators_per_organism")


tcc_object <- list(counts = abiotic_stresses$raw_counts)
threshold = 10*length(abiotic_stresses$conditions)
tcc_object$counts <- tcc_object$counts[rowSums(tcc_object$counts) > threshold,]
normalized_counts <- tcc_object$counts


fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)
topTags <- DIANE::estimateDEGs(fit, reference = "M", perturbation = "MH", p.value = 0.05, lfc = 1.5)

# adding annotations
genes <- stringr::str_split_fixed(topTags$table$genes, '\\.', 2)[,1]


tfs <- intersect(genes, regulators_per_organism$`Arabidopsis thaliana`)
targets <- genes



write.table(targets, file = "data/targets.csv", quote = F, row.names = F)
write.table(tfs, file = "data/tfs.csv", quote = F, row.names = F)
