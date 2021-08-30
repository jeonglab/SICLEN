# SICLEN : Single-cell Clustering based on effective noise reduction through ensenble similarity network


## Abstract


Single-cell sequencing provides novel means to interpret the transcriptomic profiles of individual cells. To obtain in-depth analysis of single-cell sequencing, it requires effective computational methods to accurately predict single-cell clusters because single-cell sequencing techniques only provide the transcriptomic profiles of each cell. Although the accurate estimation of the cell-to-cell similarity is an essential first step to derive reliable single-cell clustering results, it is challenging to obtain accurate similarity measurements because it highly depends on the selection of genes for similarity evaluations and the optimal set of genes for the accurate similarity estimation is typically unknown. Moreover, due to technical limitations, single-cell sequencing includes a larger number of artificial zeros and the technical noise makes it difficult to develop effective single-cell clustering algorithms. Here, we describe a novel single-cell clustering algorithm that can accurately predict single-cell clusters in large-scale single-cell sequencing by effectively reducing the zero-inflated noise and accurately estimating the cell-to-cell similarities. First, we construct an ensemble similarity network based on different similarity estimates, and reduce the artificial noise using a random walk with restart approach. Finally, starting from a larger number small but highly consistent clusters, we iteratively merge a pair of clusters with the maximum similarities until it reaches the predicted number of clusters. Extensive performance evaluation shows that the proposed single-cell clustering algorithm can yield the accurate single-cell clustering results and it can help deciphering the key messages underlying complex biological mechanisms. 


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

The sample data includes two arguments: i) gene expression counts and ii) true cell type labels.

```
> usoskin$counts
                        L128_D01  L128_H01  L128_B02  L128_C02  L128_E02  L128_A03  L128_B03  L128_D03  L128_E03  L128_C04
Xkr4                        0.00     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000
                         L128_F04  L128_G04  L128_B05  L128_D05  L128_G05  L128_D06  L128_F06  L128_H07  L128_A08  L128_H08
Xkr4                        0.000     0.000     0.000     0.000    32.027     0.000     0.000     0.000     0.000     0.000
 [ reached getOption("max.print") -- omitted 19533 rows ]
 
 > usoskin$label
  [1] "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF" 
 [20] "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF"  "NF" 
 ....
[609] "TH"  "TH"  "TH"  "TH"  "TH"  "TH"  "TH"  "TH"  "TH"  "TH"  "TH"  "TH"  "TH"  "TH" 

```

Next, run SICLEN to obtain the accurate single-cell clustering resutls. We only need to type the following command. 
```
clustering <- siclen(usoskin$counts)
```

Once the proposed single-cell clustering algorithm finished, it returns the following objects: 
```
clustering$counts      # noise reduced gene expression counts (cpm normalized)
clustering$network     # learned ensemble similarity network
clustering$clusterid   # predicted clustering labels
clustering$fgenes      # potential feature genes that are employed to construc the ensemble similarity network   
```

Finally, to verify the clustering results, we can obtain the low-dimensional visualize of single-cell clustering results through t-SNE. In this example, we will use the R packages ggplot2 and Rtnse. Please refere to the following sample code:

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



