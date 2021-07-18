# SICLEN : Single-cell Clustering based on effective noise reduction through ensenble similarity network


## Abstract


Single-cell RNA sequencing provides novel methods to interpret the transcriptomic profiles of individual cells. To obtain in-depth analysis of single-cell RNA sequencing, it requires effective computational methods to accurately predict single-cell clusters in the data because single-cell RNA sequencing techniques only provides the transcriptomic profiles of each cell. Although accurate estimation of the cell-to-cell similarity is the essential first step to derive reliable single-cell clustering, it is challenging to obtain accurate similarity measurements because it highly depends on the selection of genes for similarity evaluations and the optimal set of genes for accurate similarity estimation is typically unknown. Moreover, due to the technical limitations, single-cell RNA sequencing includes larger number of artificial zeros and these technical noise makes it difficult to develop effective single cell clustering algorithms. Here, we describe novel single-cell clustering algorithm that can accurately predict cell types in large-scale single-cell RNA sequencing by effectively reducing zero-inflated noise and accurately estimating the cell-to-cell similarities. First, we construct an ensemble similarity network based on different similarity estimates, and reduce the artificial noise using the random walk with restart approach. Finally, starting from a larger number small clusters, we iteratively merge a pair of clusters with the maximum similarities until it reaches the predicted number of clusters. Extensive performance evaluation shows that the proposed single-cell clustering algorithm can yield accurate single-cell clusters and it can help deciphering the key messages underlying complex biological mechanisms.



<p align="center">
<img src = ./vignettes/overview.png width="80%">
</p>

<p align="center">
Fig. 1 Graphical overview of SICLEN
</p>








## Installation guide 
To install SICLEN (R package), you need to install devtools and type the following command:
```
install.packages("devtools")
library(devtools)
install_github("jeonglab/SICLEN")
```



## Example 
To test the package, we will use the sample single-cell sequencing (Usoskin et al. 2015) and verify the clustering results through low-dimensional visualization



To begin with, we need to load the R package SICLEN and sample single-cell RNA sequencing data. Please type the following code:
```
library(SICLEN)
data(usoskin)
```

Next, run SICLEN to obtain the accurate single-cell clustering resutls. We only need to type the following command. 
```
clustering <- siclen(usoskin$counts)
```
Clustering resutls includes following objects: 
1) predicted clustering labels, 2) learned ensemble similarity network, 3) potential feature genes, and 4) noise reduced gene expression counts (cpm normalized)

Finally, to visualize the single-cell clustering results, we employ t-SNE. 

```
# load required packages
library(ggplot2)
library(Rtsne)

# run t-SNE
tsne <- Rtsne(t(log2(1+clustering$counts[clustering$fgenes, ])))
df.tsne <- data.frame(tsne$Y, label = factor(clustering$clusterid))

# visualization using t-SNE
ggplot(df.tsne,aes(x = X1, y = X2, col = factor(label)))+
  geom_point(size = 1, alpha = 0.7)+
  labs(x= "Dim_1", y = "Dim_2", title = 'Low-dim. visualization (t-SNE)')+
  labs(colour = "Label")


```
<p align="center">
<img src = ./vignettes/sample_tsne.png width="50%">
</p>

<p align="center">
Fig. 2 Low-dimensional visualization of the clustering resutls
</p>




## Reference
Usoskin, Dmitry, et al. "Unbiased classification of sensory neuron types by large-scale single-cell RNA sequencing." Nature Neuroscience 18.1 (2015): 145-153.





## Citation info.
Hyundoo Jeong, Sungtae Shin, and Hong Gi Yeom "Accurate single-cell clustering through effective noise reductionover ensemble similarity network" Submitted to review.



