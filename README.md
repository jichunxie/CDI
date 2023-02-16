## Clustering Deviance Index (CDI)

This repository stores the R package of the paper Clustering Deviation Index (CDI): A robust and
accurate internal measure for evaluating scRNA-seq data clustering. 

### Installation guidance

In R, run the following block of codes:

```
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("jichunxie/CDI", build_vignettes = TRUE) 


```

A tutorial of this package can be found by running the following block of codes in R:

```
library(CDI)
browseVignettes("CDI")
```
### Related paper
Clustering Deviation Index (CDI): a robust and accurate internal measure for evaluating scRNA-seq data clustering
[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02825-5](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02825-5)


