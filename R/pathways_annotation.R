#' Pathways annotation
#'
#' @useDynLib pwannot
#'
#' @import Matrix
#' @import Seurat
#' @importFrom stats p.adjust phyper
#'
#' @param object is a Seurat object after scaling and clustering
#' @param genes_list is a list of pathways which are included a pathway with a genes list, each of them contained the pathways genes in a character format (notice you should check whether genes of your organism capitalize or not)
#' @param min_length is a minimal size (number of genes) of pathways which will be analyzed
#' @param max_length is a maximal size (number of genes) of pathways which will be analyzed
#' @param p_val_border is a decimal places of p-value in cumulative matrix which define the number of success states in hypergeometric distribution
#' @param sample_number is a nubmer of random sample generations (10000 is recommended)
#'
#' @return the p-value matrix with pathways as rows and cell clusters as cols
#' @export
#'
#' @examples pathways_annotation(pbmc, genes_list, 20, 500, 4, 10000)
pathways_annotation <- function(object, genes_list, min_length, max_length, p_val_border, sample_number) {
  data <- object@scale.data
  length_of_pathways <- lapply(genes_list, function(genes) {
    length(which(rownames(data) %in% genes))
  })
  genes_list <- genes_list[length_of_pathways >= min_length & length_of_pathways <= max_length]
  length_of_pathways <- length_of_pathways[length_of_pathways >= min_length & length_of_pathways <= max_length]
  N <- length(object@data@Dimnames[[2]])
  mat <- matrix(0L, nrow = length(genes_list), ncol = ncol(data))
  pb <- progress_bar$new(
    format = "Create the reference matrix [:bar] :percent in :elapsed",
    total = length(genes_list), clear = FALSE, width= 60)
  for (i in 1:length(genes_list)){
    pb$tick()
    genes <- genes_list[[i]]
    genes_indexes <- which(rownames(data) %in% genes)
    genes_indexes <- genes_indexes - 1
    mat[i,] <- reference_matrix(data, genes_indexes)
  }
  result_pvals <- matrix(nrow = length(genes_list), ncol = length(levels(object@ident)))
  pval <- cumulative_matrix(data, unlist(length_of_pathways), mat, sample_number)
  pb <- progress_bar$new(
    format = "Hypergeometric distribution [:bar] :percent in :elapsed",
    total = nrow(pval) * length(levels(object@ident)), clear = FALSE, width= 60)
  for (l in 1:nrow(pval)) {
    K <- sum(-log10(pval[l,]) > p_val_border)
    for (j in 1:length(levels(object@ident))){
      pb$tick()
      i <- levels(object@ident)[j]
      indexes <- which(object@ident == i)
      n <- sum(object@ident == i)
      k <- sum(object@ident == i & -log10(pval[l, ]) > -log10(1/sample_number))
      result_pvals[l,j] <- phyper(k, K, N-K, n, lower.tail = FALSE)
    }
  }
  result_pvals[] <- p.adjust(result_pvals, method = 'bonferroni')
  colnames(result_pvals) <- levels(object@ident)
  rownames(result_pvals) <- names(genes_list)
  return(result_pvals)
}
' Pathways annotation
#' 
#' @useDynLib pwannot
#' 
#' @import Matrix
#' @import Seurat
#' @importFrom stats p.adjust phyper
#' 
#' @param object is a Seurat object after scaling and clustering
#' @param genes_list is a list of pathways which are included a pathway with a genes list, each of them contained the pathways genes in a character format (notice you should check whether genes of your organism capitalize or not)
#' @param min_length is a minimal size (number of genes) of pathways which will be analyzed
#' @param max_length is a maximal size (number of genes) of pathways which will be analyzed
#' @param sample_number is a nubmer of random sample generations (10000 is recommended)
#'
#' @return the p-value matrix with pathways as rows and cell clusters as cols
#' @export
#' 
#' @examples
pathways_annotation <- function(object, genes_list, min_length, max_length, sample_number) {
  data <- object@scale.data
  length_of_pathways <- lapply(genes_list, function(genes) {
    length(which(rownames(data) %in% genes))
  })
  genes_list <- genes_list[length_of_pathways >= min_length & length_of_pathways <= max_length]
  length_of_pathways <- length_of_pathways[length_of_pathways >= min_length & length_of_pathways <= max_length]
  N <- length(object@data@Dimnames[[2]])
  mat <- matrix(0L, nrow = length(genes_list), ncol = ncol(data))
  for (i in 1:length(genes_list)){
    genes <- genes_list[[i]]
    genes_indexes <- which(rownames(data) %in% genes)
    genes_indexes <- genes_indexes - 1
    mat[i,] <- reference_matrix(data, genes_indexes)
  }
  result_pvals <- matrix(nrow = length(genes_list), ncol = length(levels(object@ident)))
  pval <- cumulative_matrix(data, unlist(length_of_pathways), mat, sample_number)
  for (l in 1:nrow(pval)) {
    K <- sum(-log10(pval[l,]) > -log10(1/sample_number))
    for (j in 1:length(levels(object@ident))){
      i <- levels(object@ident)[j]
      indexes <- which(object@ident == i)
      n <- sum(object@ident == i)
      k <- sum(object@ident == i & -log10(pval[l, ]) > -log10(1/sample_number))
      result_pvals[l,j] <- phyper(k, K, N-K, n, lower.tail = FALSE)
    }
  }
  result_pvals[] <- p.adjust(result_pvals, method = 'bonferroni')
  colnames(result_pvals) <- levels(object@ident)
  rownames(result_pvals) <- names(genes_list)
  return(result_pvals)
}

