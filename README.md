# pwannot
=========

R package for automated pathway annotation in single-cell RNA-seq

# Installing the package
------------------------
You can install the package using devtools::install_github:

```{r}
devtools::install_github("mariafiruleva/pwannot")
```
# Getting started
-------------

To start working, we must have a scaled and clustered scRNA-seq dataset and the set of pathways. We will use the result from Seurat package tutorial ()

```{r}
library(pwannot)
pbmc <- readRDS('scaled_PBMC.rdata')
genes_list <- readRDS(curatedSymbol.rdata)
```
# Find upregulated pathways
---------------------------

Next run the analysis (time-consuming step). "10" is the minimal size of pathway, "500" is the maximal size of pathway, "10000" is a nubmer of random sample generations (10000 is recommended).

```{r}
annotation_results <- pathways_annotation(pbmc, genes_list, 10, 500, 10000)

```

Now we have a data frame with clusters as columns and pathway names as rows with values as adjusted p-values.

# Visualize the interested pathway
----------------------------------

For visualization the expression distribution of interested pathway (e.g., "KEGG_ALLOGRAFT_REJECTION"), we use simply function:

```{r}
plot_target_pw(pbmc, annotation_results, "KEGG_ALLOGRAFT_REJECTION")
```

![](Readme_files/figure-markdown_github/networks-1.png)



