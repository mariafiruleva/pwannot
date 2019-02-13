#' Pathways annotation
#'
#' @useDynLib pwannot
#'
#' @import Matrix
#' @import Seurat
#' @import progress
#' @importFrom stats p.adjust phyper
#'
#' @param object is a Seurat object after scaling and clustering
#' @param genes_list is a list of pathways which are included a pathway with a genes list, each of them contained the pathways genes in a character format (notice you should check whether genes of your organism capitalize or not)
#' @param min_length is a minimal size (number of genes) of pathways which will be analyzed
#' @param max_length is a maximal size (number of genes) of pathways which will be analyzed
#' @param significance_level is a decimal places of p-value in cumulative matrix which define the number of success states in hypergeometric distribution
#' @param sample_number is a nubmer of random sample generations (10000 is recommended)
#'
#' @return the p-value matrix with pathways as rows and cell clusters as cols
#' @export
#'
#' @examples pathways_annotation(pbmc, genes_list, 20, 500, 0.05, 10000)
pathways_annotation <- function(object, genes_list, min_length = 10, max_length = 500,
                                significance_level = 0.05, sample_number = 1000) {
 
  data <- GetAssayData(object, slot = "scale.data")
  N <- length(GetAssayData(object, slot = "data")@Dimnames[[2]])
 
  length_of_pathways <- lapply(genes_list, function(genes) {
    sum(genes %in% rownames(data))
  })
 
  genes_list <- genes_list[length_of_pathways >= min_length & length_of_pathways <= max_length]
  length_of_pathways <- length_of_pathways[length_of_pathways >= min_length & length_of_pathways <= max_length]
 
  message(sprintf("%d pathways left after filtering", length(genes_list)))
 
  mat <- matrix(0L, nrow = length(genes_list), ncol = ncol(data))
 
  pb <- progress_bar$new(
    format = "Creating the reference matrix [:bar] :percent in :elapsed",
    total = length(genes_list), clear = FALSE, width=60, force=TRUE)
 
  for (i in 1:length(genes_list)){
    pb$tick()
    genes <- genes_list[[i]]
    genes_indexes <- which(rownames(data) %in% genes)
    genes_indexes <- genes_indexes - 1
    mat[i,] <- reference_matrix(data, genes_indexes)
  }
  result_pvals <- matrix(nrow = length(genes_list), ncol = length(levels(Idents(object))))
  message("Sampling starts")
  pval <- cumulative_matrix(data, unlist(length_of_pathways), mat, sample_number)
  pb <- progress_bar$new(
    format = "Getting hypergeometric p-values [:bar] :percent in :elapsed",
    total = nrow(pval) * length(levels(Idents(object))), clear = FALSE, width= 60)
  for (l in 1:nrow(pval)) {
    K <- sum(pval[l, ] < significance_level)
    for (j in 1:length(levels(Idents(object)))){
      pb$tick()
      i <- levels(Idents(object))[j]
      indexes <- which(Idents(object) == i)
      n <- sum(Idents(object) == i)
      k <- sum(Idents(object) == i & pval[l, ] < significance_level)
      result_pvals[l,j] <- phyper(k, K, N-K, n, lower.tail = FALSE)
    }
  }
  result_pvals[] <- p.adjust(result_pvals, method = 'bonferroni')
  colnames(result_pvals) <- levels(Idents(object))
  rownames(result_pvals) <- names(genes_list)
  return(result_pvals)
}

