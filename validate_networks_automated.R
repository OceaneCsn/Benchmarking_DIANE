library(stringr)
library(igraph)
# ______ Validation of inferred networks from DIANE's thresholding procedure _________ #


# ___ Generates expression data to use for benchmark ______ #

library(DIANE)
data("abiotic_stresses")
data("gene_annotations")
data("regulators_per_organism")


tcc_object <- list(counts = abiotic_stresses$raw_counts)
threshold = 10*length(abiotic_stresses$conditions)
tcc_object$counts <- tcc_object$counts[rowSums(tcc_object$counts) > threshold,]
normalized_counts <- tcc_object$counts


fit <- DIANE::estimateDispersion(tcc = tcc_object, conditions = abiotic_stresses$conditions)



# ___ Get a predicted network with given settings ______ #

get_network <- function(nEdges = NULL, lfc = 2, nTrees = 1000, 
                        nShuffle = 1000,
                        nCores = 20, density = 0.03, fdr = 0.05){
  
  topTags <- DIANE::estimateDEGs(fit, reference = "C", perturbation = "H", p.value = 0.05, lfc = lfc)
  genes <- get_locus(topTags$table$genes)
  regressors <- intersect(genes, regulators_per_organism[["Arabidopsis thaliana"]])
  aggregated_data <- aggregate_splice_variants(data = normalized_counts)
  
  grouping <- DIANE::group_regressors(aggregated_data, genes, regressors)
  
  grouped_counts <- grouping$counts
  grouped_targets <- grouping$grouped_genes
  grouped_regressors <- grouping$grouped_regressors
  
  
  mat <- DIANE::network_inference(grouped_counts, 
                                  conds = abiotic_stresses$conditions, 
                                  targets = grouped_targets, 
                                  regressors = grouped_regressors, 
                                  importance_metric = "MSEincrease_oob",
                                  nCores = nCores, verbose = FALSE, nTrees = nTrees)
  
  nGenes = length(grouped_targets)
  nRegulators = length(grouped_regressors)
  
  
  if(!is.null(nEdges)){
    network <- DIANE::network_thresholding(mat, n_edges = nEdges)
    d <- network_data(network, 
                      regulators_per_organism[["Arabidopsis thaliana"]], 
                      gene_annotations$`Arabidopsis thaliana`)
  }
  else{
    res <- DIANE::test_edges(mat, normalized_counts = grouped_counts, density = density,
                             nGenes = nGenes, 
                             nRegulators = nRegulators, 
                             nTrees = nTrees, nCores = nCores)
    
    network <- DIANE::network_from_tests(res$links, fdr = fdr)
    d <- network_data(network, 
                      regulators_per_organism[["Arabidopsis thaliana"]], 
                      gene_annotations$`Arabidopsis thaliana`)
  }
  return(list(edges = d$edges, nGenes = nGenes, nRegulators = nRegulators))
}



 #net <- get_network(270, nTrees = 2000)

#' Deals with grouped nodes of a network
#' 
#' from mean_ATGSJHKS-ATGODLKSL to ATGDHJKHK 
#' becomes
#' from ATGSJHKS to ATGDHJKHK
#' from ATGODLKSL to ATGDHJKHK
#' 
#' @param net network containing grouped regulators,
#' in the form of a dataframe of edges
#'
#' @return dataframe of edges with duplicated edges for grouped regulators
flatten_edges <- function(net){
  distinct <- net[!str_detect(net$from, 'mean_') & !str_detect(net$to, 'mean_'),c("from", "to")]
  grouped <- net[str_detect(net$from, 'mean_') | str_detect(net$to, 'mean_'),]
  
  if(nrow(grouped) > 0){
    for(i in 1:nrow(grouped)){
      tf <- grouped[i,"from"]
      targ <- grouped[i,"to"]
      if(str_detect(tf, "mean_")){
        tfs <- unlist(strsplit(stringr::str_remove(tf, "mean_"), '-'))
      }
      else
        tfs <- tf
      if(str_detect(targ, "mean_")){
        targs <- unlist(strsplit(stringr::str_remove(targ, "mean_"), '-'))
      }
      else 
        targs <- targ
      for(tfi in c(tfs)){
        for(targi in c(targs)){
          distinct <- rbind.data.frame(distinct, c(tfi, targi))
        }
      }
    }
  }
  
  return(distinct)
}




# compute fraction of validated edges
validate_network <- function(net, no_dap = FALSE){
  
  # loading edges validated by connecTF database
  connecTF <- read.csv("data/connectf_CH.csv")
  connecTF <- connecTF[,c("gene_id", "TARGET", "TECHNOLOGY.METHOD", "EXPERIMENT_TYPE")]
  if(no_dap)
    connecTF <- connecTF[connecTF$TECHNOLOGY.METHOD != "DAPSeq",]
  validated <- igraph::graph_from_data_frame(connecTF, directed = TRUE, vertices = NULL)
  table(connecTF$TECHNOLOGY.METHOD)

  # ungroup grouped regulators to validate individual interactions
  flat <- flatten_edges(net)
  pred <- graph_from_data_frame(flat, directed = TRUE, vertices = NULL)
  
  # Number of edges for which we have validation data
  tfs_val <- intersect(flat$from, connecTF$gene_id)
  validable_edges <- nrow(flat[flat$from %in% tfs_val,])
  inter <- intersection(pred, validated, keep.all.vertices = TRUE)
  validated_edges <- gsize(inter)
  validation_rate <- validated_edges/validable_edges
  print(paste(validated_edges, "edges were validated over the", validable_edges, "edges that had validation information."))
  print(paste("Validation rate :", round(validation_rate*100, 2), "%"))
  return(validation_rate)
}

# N <- 10
# nCores <- 32
# nTrees <- 2000
# 
# lfcs <- c(2, 1.5)
# densities <- c(0.03, 0.01)
# fdrs <- c(0.05, 0.01)
# strategies <- c("testing", "before_testing", "same_edges")
# 
# results <- setNames(expand.grid(strategies, densities, fdrs, lfcs), c("Strategy", "density", "fdr", "lfc"))
# results <- results[rep(seq_len(nrow(results)), each = N), ]
# results$precision <- NA



# ___ Starts the benchmark ___ #

benchmark <- function(nTrees = 1000, nCores = 32, N = 15,
                      nShuffle = 1000, no_dap = FALSE,
                      outfile = "benchmark_1000Trees_CvsH_withDap.csv"){
  lfcs <- c(2, 1.5)
  densities <- c(0.03, 0.01)
  
  close(file(outfile, open="w"))
  to_store <- c("lfc", "density", "fdr", "strategy", "replicate", "precision")
  write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
  
  for(lfc in lfcs){
    for(density in densities){
      for(i in 1:N){
        
        # testing procedure
        tmp <- get_network(lfc = lfc, nTrees = nTrees, nCores = nCores, 
                           density = density, fdr = 0.05, nShuffle = nShuffle)
        net <- tmp$edges
        nGenes <- tmp$nGenes
        nRegulators <- tmp$nRegulators
        n_before <- DIANE::get_nEdges(density = density, 
                                      nGenes = nGenes, 
                                      nRegulators = nRegulators)
        n_after <- nrow(net)
        net_0.01 <- net[net$fdr <= 0.01,]
        n_after_0.01 <- nrow(net_0.01)
        
        
        ## fdr 0.05
        precision <- validate_network(net, no_dap =no_dap)
        to_store <- c(lfc, density, 0.05, "testing", i, precision)
        write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
        
        
        
        # network before testing
        tmp <- get_network(nEdges = n_before, lfc = lfc, nTrees = nTrees, 
                           nCores = nCores, density = density, fdr = 0.05)
        
        precision <- validate_network(tmp$edges, no_dap =no_dap)
        to_store <- c(lfc, density, 0.05, "before_testing", i, precision)
        write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
        
        # network same edges
        tmp <- get_network(nEdges = n_after, lfc = lfc, nTrees = nTrees, 
                           nCores = nCores, density = density, fdr = 0.05)
        
        precision <- validate_network(tmp$edges, no_dap =no_dap)
        to_store <- c(lfc, density, 0.05, "same_edges", i, precision)
        write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
        
        ## fdr 0.01
        
        precision <- validate_network(net_0.01, no_dap =no_dap)
        to_store <- c(lfc, density, 0.01, "testing", i, precision)
        write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
        
        # network before testing
        tmp <- get_network(nEdges = n_before, lfc = lfc, nTrees = nTrees, 
                           nCores = nCores, density = density, fdr = 0.05)
        
        precision <- validate_network(tmp$edges, no_dap =no_dap)
        to_store <- c(lfc, density, 0.01, "before_testing", i, precision)
        write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
        
        # network same edges
        tmp <- get_network(nEdges = n_after_0.01, lfc = lfc, nTrees = nTrees, 
                           nCores = nCores, density = density, fdr = 0.05)
        
        precision <- validate_network(tmp$edges, no_dap =no_dap)
        to_store <- c(lfc, density, 0.01, "same_edges", i, precision)
        write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      }
    }
  }
}


benchmark(nTrees = 1000, nCores = 40, N = 8, nShuffle = 1000, no_dap = FALSE)


data <- read.csv("benchmark_1000Trees_CvsH_withDap.csv")


library(ggplot2)
library(ggpubr)


ggplot(data, aes(color = strategy, fill = strategy, x = strategy, y = precision)) + 
  geom_boxplot(size = 0.5, alpha = 0.3) + geom_jitter(size = 0.5)+ ylim(0.1,0.4) +
  facet_wrap(~lfc + density + factor(fdr), nrow = 2) + scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") + stat_compare_means(
    aes(x = strategy, y = precision),
    comparisons = list(c("testing", "same_edges"), 
                       c("testing", "before_testing")), method = "wilcox.test", paired = FALSE) + 
  ggtitle("Precision on connecTF, C vs H genes, nTrees = 1000, nShuffle = 1000, N = 15, with dap seq")

 

# idees : augmenter nShuffle, changer liste de genes, pourquoi lfc 1.5, 3, 0.01?, valider avec Chip

# validate_network(density = 0.03, lfc = 1.5)
# 
# res <- 
# mean(res)
# sd(res)
# 
# res <- sapply(1:N, validate_network, n_edges = 1418)
# mean(res)
# sd(res)
# 
# validate_network(density = 0.03, lfc = 1.5, fdr = 0.01)
# res <- sapply(1:N, validate_network, n_edges = 239)
# mean(res)
# sd(res)
# validate_network(density = 0.03, lfc = 1.5, fdr = 0.01, strat = "same_n_edges")
# 
# validate_network(density = 0.01, lfc = 1.5)
# validate_network(density = 0.01, lfc = 1.5, strat = "no_tests")
# validate_network(density = 0.01, lfc = 1.5, strat = "same_n_edges")
# 
# validate_network(density = 0.01, lfc = 1.5, fdr = 0.01)
# validate_network(density = 0.01, lfc = 1.5, fdr = 0.01, strat = "same_n_edges")
# 
# 
# 
# 
# validate_network(density = 0.03, lfc = 2)
# validate_network(density = 0.03, lfc = 2, strat = "no_tests")
# validate_network(density = 0.03, lfc = 2, strat = "same_n_edges")
# 
# 
# validate_network(density = 0.03, lfc = 2, fdr = 0.01)
# validate_network(density = 0.03, lfc = 2, fdr = 0.01, strat = "same_n_edges")
# 
# validate_network(density = 0.01, lfc = 2)
# validate_network(density = 0.01, lfc = 2, strat = "no_tests")
# validate_network(density = 0.01, lfc = 2, strat = "same_n_edges")
# 
# validate_network(density = 0.01, lfc = 2, fdr = 0.01)
# validate_network(density = 0.01, lfc = 2, fdr = 0.01, strat = "same_n_edges")
# 
# read.csv("data/network_edges_d_0.03_lfc_1.5_no_tests.csv")


