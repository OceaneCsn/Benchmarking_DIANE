library(stringr)
library(igraph)

# ______ Validation of inferred networks from DIANE's thresholding procedure _________ #

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
  return(distinct)
}

# loading edges validated by connecTF database
connecTF <- read.csv("data/connectf.csv")
connecTF <- connecTF[,c("gene_id", "TARGET", "TECHNOLOGY.METHOD", "EXPERIMENT_TYPE")]
validated <- igraph::graph_from_data_frame(connecTF, directed = TRUE, vertices = NULL)
table(connecTF$TECHNOLOGY.METHOD)


# compute fraction of validated edges
validate_network <- function(fdr = 0.05, density = 0.03, lfc = 1.5, strat = "testing"){
  if(strat == "testing"){
    to_read <- paste0("data/network_edges_d", density, "_lfc_", lfc, ".csv")
    net <- read.csv(to_read)
    net <- net[net$fdr < fdr,]
    print(nrow(net))
  }
    
  if(strat == "no_tests"){
    to_read <- paste0("data/network_edges_d", density, "_lfc_", lfc, "_no_tests.csv")
    net <- read.csv(to_read)
  }
    
  if(strat == "same_n_edges"){
    to_read <- paste0("data/network_edges_d", density, "_lfc_", lfc, "_same_n_edges.csv")
    if(fdr == 0.01)
      to_read <- paste0("data/network_edges_d", density, "_lfc_", lfc, "_same_n_edges_fdr0.01.csv")
    net <- read.csv(to_read)
  }
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
}


validate_network(density = 0.03, lfc = 1.5)
validate_network(density = 0.03, lfc = 1.5, strat = "no_tests")
validate_network(density = 0.03, lfc = 1.5, strat = "same_n_edges")


validate_network(density = 0.03, lfc = 1.5, fdr = 0.01)
validate_network(density = 0.03, lfc = 1.5, fdr = 0.01, strat = "same_n_edges")

validate_network(density = 0.01, lfc = 1.5)
validate_network(density = 0.01, lfc = 1.5, strat = "no_tests")
validate_network(density = 0.01, lfc = 1.5, strat = "same_n_edges")

validate_network(density = 0.01, lfc = 1.5, fdr = 0.01)
validate_network(density = 0.01, lfc = 1.5, fdr = 0.01, strat = "same_n_edges")




validate_network(density = 0.03, lfc = 2)
validate_network(density = 0.03, lfc = 2, strat = "no_tests")
validate_network(density = 0.03, lfc = 2, strat = "same_n_edges")


validate_network(density = 0.03, lfc = 2, fdr = 0.01)
validate_network(density = 0.03, lfc = 2, fdr = 0.01, strat = "same_n_edges")

validate_network(density = 0.01, lfc = 2)
validate_network(density = 0.01, lfc = 2, strat = "no_tests")
validate_network(density = 0.01, lfc = 2, strat = "same_n_edges")

validate_network(density = 0.01, lfc = 2, fdr = 0.01)
validate_network(density = 0.01, lfc = 2, fdr = 0.01, strat = "same_n_edges")

read.csv("data/network_edges_d_0.03_lfc_1.5_no_tests.csv")


