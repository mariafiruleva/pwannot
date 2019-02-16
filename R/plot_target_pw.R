#' Plot a target pathway distribution
#'
#' @useDynLib pwannot
#'
#' @import ggplot2
#'
#' @param object is a Seurat object after scaling and clustering
#' @param annotation_result is your result of pathways_annotation function
#' @param pw_name is a name of interested pathway
#'
#' @return the plot of a pathway distribution
#' @export
#'
#' @examples plot_target__pw(pbmc, annotation_result, "KEGG_ALLOGRAFT_REJECTION")
plot_target_pw <- function(object, gene_set, reduction="tsne") {
  data <- GetAssayData(object, slot = "scale.data")
  tsne <- object@reductions[[reduction]]@cell.embeddings
  tsne <- as.data.frame(tsne)
  colnames(tsne) <- paste0(reduction, 1:ncol(tsne))
 
  genes_indexes <- which(rownames(data) %in% gene_set[[1]])
  expression_signaling <- colMeans(data[genes_indexes,])
  tsne$expression <- expression_signaling
 
  ggplot(tsne, aes_string(x=paste0(reduction, 1),
                          y=paste0(reduction, 2),
                          color="expression"))+
    geom_point()+
    scale_color_gradient2(low='blue', high='red', mid = 'grey')+
    theme_classic(base_size=8) +
    theme(aspect.ratio = 1)
}