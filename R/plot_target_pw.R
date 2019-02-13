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
plot_target_pw <- function(object, annotation_result, pw_name) {
  data <- object@scale.data
  tsne <- object@dr$tsne@cell.embeddings
  tsne <- as.data.frame(tsne)
  pathway_name <- which(grepl(pw_name, rownames(annotation_result), ignore.case = TRUE))
  pathway_name <- as.character(as.name(row.names(annotation_result)[pathway_name]))
  pathway_genes <- curatedSymbol[pathway_name]
  genes_indexes <- which(rownames(object@scale.data) %in% pathway_genes[[1]])
  expression_signaling <- colMeans(object@scale.data[genes_indexes,])
  tsne$expression <- expression_signaling
  ggplot(tsne, aes(x=tSNE_1, y=tSNE_2, color=expression))+
    geom_point()+
    scale_color_gradient2(low='blue', high='red', mid = 'grey')+
    theme_classic(base_size=8)+
    theme(aspect.ratio = 1)+
    ggtitle(pathway_name)
}

