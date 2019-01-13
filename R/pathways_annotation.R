#' Title
#'
#' @param scaled_object is a scaled Seurat object (seurat_object@scale.data))
#' @param genes_list is a list of pathways which are included a pathway with a genes list, each of them contained the pathways genes in a character format (notice you should check whether genes of your organism capitalize or not)
#' @param sample_number is a nubmer of random sample generations (10000 is recommended)
#'
#' @return the p-value matrix with pathways as rows and cell clusters as cols
#' @export
#'
#' @examples
#' pathways_annotation(seurat_object, genes_list, 10000)
pathways_annotation <- function(object, genes_list, sample_number) {
  data <- object@scale.data
  length_of_pathways <- lapply(genes_list, function(genes) {
    length(which(rownames(data) %in% genes))
  })
  genes_list <- genes_list[length_of_pathways >= 10]
  length_of_pathways <- length_of_pathways[length_of_pathways >= 10]
  N <- length(object@data@Dimnames[[2]])
  mat <- matrix(0L, nrow = length(genes_list), ncol = ncol(data))
  for (i in 1:length(genes_list)){
    genes <- genes_list[[i]]
    genes_indexes <- which(rownames(data) %in% genes)
    genes_indexes <- genes_indexes - 1
    mat[i,] <- reference_matrix(data, genes_indexes)
  }
  result_pvals <- matrix(nrow = length(genes_list), ncol = length(levels(object@ident)))
  expression_signaling <- colMeans(data[genes_indexes,])
  genes_indexes <- genes_indexes - 1
  pval <- cumulative_matrix(data, unlist(length_of_pathways), mat, sample_number)
  for (l in 1:nrow(pval)) {
    K <- sum(-log10(pval[l,]) > 2)
    for (j in 1:length(levels(object@ident))){
      i <- levels(object@ident)[j]
      indexes <- which(object@ident == i)
      n <- sum(object@ident == i)
      k <- sum(object@ident == i & -log10(pval[l, ]) > 2)
      result_pvals[l,j] <- phyper(k, K, N-K, n, lower.tail = FALSE)
    }
  }
  result_pvals[] <- p.adjust(result_pvals, method = 'bonferroni')
  colnames(result_pvals) <- levels(object@ident)
  rownames(result_pvals) <- names(genes_list)
  return(result_pvals)
}