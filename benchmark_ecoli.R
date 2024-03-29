library(DIANE)
library(tidyverse)
library(tictoc)
library(igraph)
library(tigress)

data("ecoli")

nCores = 40
nTrees = 1000

 

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
  distinct <- net[!str_detect(net$from, 'mean_') & !str_detect(net$to, 'mean_'), c("from", "to")]
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



get_mat <- function(i, nCores = 40, nTrees = 1000){
  
  data <- ecoli$exp
  #d <-ecoli$reg
  data <- data[,str_detect(colnames(data), 'T0|T24|T48|T12|T36|T60')]
  colnames(data) <- str_replace(colnames(data), "_N", ".N")
  
  regulators <- rownames(data)[unique(ecoli$reg$tf)]

  genes <- rownames(data)
  
  grouping <- DIANE::group_regressors(data.frame(data), genes, regulators)
  
  grouped_counts <- grouping$counts
  grouped_targets <- grouping$grouped_genes
  grouped_regressors <- grouping$grouped_regressors
  
  conds <- str_split_fixed(colnames(grouped_counts), '_',2)[,1]
  
  print("Running GENIE3")
  print(dim(grouped_counts))
  
  tic("GENIE3 inference")
  mat <- DIANE::network_inference(grouped_counts, 
                                  conds = conds,
                                  targets = grouped_targets, 
                                  regressors = grouped_regressors, 
                                  importance_metric = "MSEincrease_oob",
                                  nCores = nCores, verbose = FALSE, nTrees = nTrees)
  toc()
  return(list(mat = mat, grouping = grouping))
}


get_network_from_nEdges <- function(mat, nEdges){
  
  network <- DIANE::network_thresholding(mat, n_edges = nEdges)
  d <- network_data(network, 
                    regulators_per_organism[["Arabidopsis thaliana"]], 
                    gene_annotations$`Arabidopsis thaliana`)
  return(d$edges)
}


get_network_before_testing <- function(mat, grouping, density = 0.01){
  
  grouped_targets <- grouping$grouped_genes
  grouped_regressors <- grouping$grouped_regressors
  
  nGenes = length(grouped_targets)
  nRegulators = length(grouped_regressors)
  
  nEdges <- get_nEdges(density, nGenes, nRegulators)
  
  return(get_network_from_nEdges(mat, nEdges))
  
}

get_network_testing <- function(mat, grouping, density,
                                nTrees = 1000, nCores = 32){
  
  grouped_counts <- grouping$counts
  grouped_targets <- grouping$grouped_genes
  grouped_regressors <- grouping$grouped_regressors
  
  nGenes = length(grouped_targets)
  nRegulators = length(grouped_regressors)
  
  res <- DIANE::test_edges(mat, normalized_counts = grouped_counts, 
                           density = density,
                           nGenes = nGenes, 
                           nRegulators = nRegulators, 
                           nTrees = nTrees, nCores = nCores)
  
  
  d_0.1 <- network_data(DIANE::network_from_tests(res$links, fdr = 0.1), 
                        regulators_per_organism[["Arabidopsis thaliana"]], 
                        gene_annotations$`Arabidopsis thaliana`)
  
  d_0.05 <- network_data(DIANE::network_from_tests(res$links, fdr = 0.05), 
                         regulators_per_organism[["Arabidopsis thaliana"]], 
                         gene_annotations$`Arabidopsis thaliana`)
  
  d_0.01 <- network_data(DIANE::network_from_tests(res$links, fdr = 0.01), 
                         regulators_per_organism[["Arabidopsis thaliana"]], 
                         gene_annotations$`Arabidopsis thaliana`)
  
  d_0.005 <- network_data(DIANE::network_from_tests(res$links, fdr = 0.005), 
                         regulators_per_organism[["Arabidopsis thaliana"]], 
                         gene_annotations$`Arabidopsis thaliana`)
  
  return(list(edges_0.1 = d_0.1$edges,
              edges_0.05 = d_0.05$edges,
              edges_0.01 = d_0.01$edges,
              edges_0.005 = d_0.005$edges))
  
}



validate_network <- function(net){
  # loading edges validated
  gold <- ecoli$reg
  
  data <- ecoli$exp
  gold$tf <- rownames(data)[gold$tf]
  gold$target <- rownames(data)[gold$target]
  
  validated <- graph_from_data_frame(gold, directed = TRUE, vertices = NULL)
  
  # ungroup grouped regulators to validate individual interactions
  flat <- flatten_edges(net)
  pred <- graph_from_data_frame(flat, directed = TRUE, vertices = NULL)
  
  # Number of edges for which we have validation data
  tfs_val <- intersect(flat$from, gold$tf)
  validable_edges <- nrow(flat[flat$from %in% tfs_val,])
  
  inter <- intersection(pred, validated, keep.all.vertices = FALSE)
  validated_edges <- gsize(inter)
  validation_rate <- validated_edges/validable_edges
  
  recall <- validated_edges/nrow(gold)
  
  print(paste(validated_edges, "edges were validated over the", validable_edges, "edges that had validation information."))
  print(paste("Validation rate :", round(validation_rate*100, 2), "%"))
  return(c(validation_rate, recall, nrow(net)))
}

  
benchmark <- function(i, nTrees = 1000, nCores = 32, N = 15, nShuffle = 1000,
                      outfile = paste0("benchmark_ecoli_", nTrees, "Trees_",N,"N.csv")){
  
  densities <- c(0.005, 0.0025)
  
  close(file(outfile, open="w"))
  to_store <- c("density", "fdr", "Strategy", "replicate", "precision", "recall", "N_edges")
  write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
  
  for(density in densities){
    for(j in 1:N){
      
      # get GENIE3 ranking
      genie3 <- get_mat(i, nCores = nCores, nTrees = nTrees)
      mat <- genie3$mat
      grouping <- genie3$grouping
      
      # testing procedure
      edges <- get_network_testing(mat = mat, grouping = grouping, 
                                   density = density, nTrees = nTrees,
                                   nCores = nCores)
      
      # fdr 0.1
      to_store <- c(density, 0.1, "testing", j, validate_network(edges$edges_0.1))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      edges_before <- get_network_before_testing(mat = mat, grouping = grouping, density = density)
      
      to_store <- c(density, 0.1, "before testing", j, validate_network(edges_before))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      edges_after <- get_network_from_nEdges(mat = mat, nEdges = nrow(edges$edges_0.1))
      to_store <- c(density, 0.1, "same n edges", j, validate_network(edges_after))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      # fdr 0.05
      to_store <- c(density, 0.05, "testing", j, validate_network(edges$edges_0.05))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      edges_before <- get_network_before_testing(mat = mat, grouping = grouping, density = density)
      
      to_store <- c(density, 0.05, "before testing", j, validate_network(edges_before))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      edges_after <- get_network_from_nEdges(mat = mat, nEdges = nrow(edges$edges_0.05))
      to_store <- c(density, 0.05, "same n edges", j, validate_network(edges_after))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      # fdr 0.01
      to_store <- c(density, 0.01, "testing", j, validate_network(edges$edges_0.01))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      edges_before <- get_network_before_testing(mat = mat, grouping = grouping, density = density)
      
      to_store <- c(density, 0.01, "before testing", j, validate_network(edges_before))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      edges_after <- get_network_from_nEdges(mat = mat, nEdges = nrow(edges$edges_0.01))
      to_store <- c(density, 0.01, "same n edges", j, validate_network(edges_after))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      # fdr 0.005
      to_store <- c(density, 0.005, "testing", j, validate_network(edges$edges_0.005))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      edges_before <- get_network_before_testing(mat = mat, grouping = grouping, density = density)
      
      to_store <- c(density, 0.005, "before testing", j, validate_network(edges_before))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
      edges_after <- get_network_from_nEdges(mat = mat, nEdges = nrow(edges$edges_0.005))
      to_store <- c(density, 0.005, "same n edges", j, validate_network(edges_after))
      write(paste(to_store, collapse = ','), file=outfile, append=TRUE)
      
    }
  }
}


benchmark(nTrees = 1000, nCores = 24, N = 10)