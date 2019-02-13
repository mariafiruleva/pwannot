# pwannot

R package for automated pathway annotation in single-cell RNA-seq

# Installing the package

You can install the package using devtools::install_github:

```{r}
devtools::install_github("mariafiruleva/pwannot")
```
# Getting started

To start working, we must have a scaled and clustered scRNA-seq dataset and the set of pathways. 

For demonstration of package functionality, we will use a dataset of 996 Peripheral Blood Mononuclear Cells (PBMCs) made publically available by 10X Genomics ([data](http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.tar.gz)) and prepared by Seurat standard workflow v3.0 [Seurat workflow]((https://satijalab.org/seurat/essential_commands.html)).

```{r}
library(pwannot)
pbmc <- readRDS('pbmc_after_seurat.rds')
genes_list <- readRDS(KEGG_pathways.rds)
```
# Find upregulated pathways

Next run the analysis (time-consuming step). "10" is the minimal size of pathway, "500" is the maximal size of pathway, "0.05" is a significance level which define the number of success states in hypergeometric distribution, "10000" is a nubmer of random sample generations (10000 is recommended).

```{r}
annotation_results <- pathways_annotation(pbmc, genes_list, 20, 500, 0.05, 10000)
```

Now we have a data frame with cluster as columns, pathways as rows with adjusted p-values.

# Visualize

For vizualisation the expression of target pathway (e.g., ""LEE_DIFFERENTIATING_T_LYMPHOCYTE"") distribution we use simply function:

```{r}
plot_target_pw(pbmc, genes_list[""LEE_DIFFERENTIATING_T_LYMPHOCYTE"], "tsne")+
ggtitle("LEE_DIFFERENTIATING_T_LYMPHOCYTE")
```


![](https://github.com/mariafiruleva/pwannot/blob/master/Readme_files/figure-markdown_github/pwannot_illustation.png)
