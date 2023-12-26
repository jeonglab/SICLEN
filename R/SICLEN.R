#' SICLEN: Single-cell Clustering based on effective noise reduction through ensenble similarity network
#' @param data (M x N) dimensional matrix for gene expressions
#' M: The number of genes
#' N: The number of cells
#' @param nens The number of similarity estimations
#' @param knum The number of nearest neighbors for KNN networks
#' @param alpha Restart probability of the random walker
#' @param npc The number of principal components
#'
#' @author Hyundoo Jeong <hdj@inu.ac.kr>
#'
#' @references
#' Accurate single-cell clustering through effective noise reduction over ensemble similarity network
#'
#' @examples
#' library(SICLEN)
#' data(scdata)
#'
#' logsc <- log2(1+cpm(scdata$counts))
#' cl_res <- SICLEN(data = logsc)
#' @export
siclen <- function(data, nens = 20, knum = 30, alpha = 0.7, npc = 10){

  print("Start: denoising")
  denoised_data <- siclen_denoise(data, nens = nens, knum = knum, alpha = alpha, npc = npc)
  print("Complete: denoising")

  log_scdata <- log2(1+cpm(denoised_data$corrected_cpm))
  myvar_genes <- IdentifyFeatures(log_scdata, pFeat = 0.01)
  var_genes <- myvar_genes$VarGenes


  print("Start: estimating the number of true clusters")
  best_nc <- est_cls(denoised_data$corrected_cpm)
  print("Complete: estimating the number of true clusters")

  print("Start: clustering")
  cl_res <- siclen_clust(denoised_data$corrected_cpm, best_nc)
  print("Complete: clustering")

  output <- list()
  output$counts <- denoised_data$corrected_cpm
  output$network <- denoised_data$eknnnet
  output$clusterid <- cl_res
  output$fgenes <- var_genes

  print("Done!!")
  return(output)

}


#' SICLEN: Single-cell Clustering based on effective noise reduction through ensenble similarity network
#' @param data (M x N) dimensional matrix for gene expressions
#' M: The number of genes
#' N: The number of cells
#' @param best_nc The number of estimated clustes
#' @param Knum Initial number of clusters
#'
#' @author Hyundoo Jeong <hdj@inu.ac.kr>
#'
#' @references
#' Accurate single-cell clustering through effective noise reduction over ensemble similarity network
#' @export
#'
siclen_clust <- function(data, best_nc, Knum = 20){
  log_scdata <- log2(1+cpm(data))
  myvar_genes <- IdentifyFeatures(log_scdata, pFeat = 0.01)
  var_genes <- myvar_genes$VarGenes

  pca <- prcomp_irlba(t(log_scdata[var_genes,]), n = 10)

  km <- kmeans(pca$x, centers = Knum)


  sind <- which(table(km$cluster) == 1)
  if(length(sind) != 0){
    print("removing singleton nodes")
    clid <- setdiff(names(table(km$cluster)), names(sind))
    for(ss in names(sind)){
      single <- which(km$cluster == ss)
      mscore <- matrix(0, ncol = length(clid), nrow=1)
      ii <- 1
      for(cc in clid){
        cid <- which(km$cluster == cc)
        mscore[ii] <- cor(log_scdata[var_genes, single], rowMeans2(log_scdata[var_genes, cid]))
        ii <- ii+1
      }
      km$cluster[single] <- clid[which(mscore == max(mscore))]
    }

  }



  memid <- km$cluster
  while(length(table(memid)) > best_nc){

    cormat <- matrix(0, nrow = Knum,  ncol = Knum)
    for (ix in 1: (Knum-1)){
      for(iy in (ix+1):Knum){
        idx <- which(memid == ix)
        idy <- which(memid == iy)
        cormat[ix, iy] <- cor(rowMeans2(log_scdata[var_genes, idx]), rowMeans2(log_scdata[var_genes, idy]))
      }
    }

    mind <- which(max(cormat) == cormat, arr.ind = T)



    idx <- which(memid == mind[2])
    memid[idx] <- mind[1]

    Knum <- Knum - 1
    memid_new <- memid
    newid <- 1
    for(iz in names(table(memid))){
      memid_new[which(memid == iz)] <- newid
      newid <- newid + 1
    }
    memid <- memid_new

  }



  return(memid)

}


#' SICLEN: Single-cell Clustering based on effective noise reduction through ensenble similarity network
#' @description
#' This function estimates the number of clusters in the single-cell RNA sequencing data
#'
#' @param data (M x N) dimensional matrix for gene expressions
#' M: The number of genes
#' N: The number of cells
#'
#' @return
#' Estimated number of clusters in the single-cell RNA sequencing data
#'
#' @author Hyundoo Jeong <hdj@inu.ac.kr>
#'
#' @references
#' Accurate single-cell clustering through effective noise reduction over ensemble similarity network
#'

#' @export
est_cls <- function(data){


  corrected_log <- log2(1+cpm(data))

  myvar_genes <- IdentifyFeatures(corrected_log, pFeat = 0.01)
  var_genes <- myvar_genes$VarGenes
  my_pca_log <- prcomp_irlba(t(corrected_log[var_genes,]), n = 10)
  PCs <- t(my_pca_log$x)

  res <- NbClust(t(PCs), distance = "euclidean",
                 method = "ward.D", index = "rubin")

  return(res$Best.nc[1])
}








#' SICLEN: Single-cell Clustering based on effective noise reduction through ensenble similarity network
#' @description
#' This function reduces the noise in the single-cell RNA sequencing through a random walk restart over the ensenble similarity network  the number of clusters in the single-cell RNA sequencing data
#'
#' @param data (M x N) dimensional matrix for gene expressions
#' M: The number of genes
#' N: The number of cells
#' @param nens The number of similarity estimations
#' @param knum The number of nearest neighbors for KNN networks
#' @param alpha Restart probability of the random walker
#' @param npc The number of principal components
#'
#' @return
#' Noise reduced single-cell RNA sequencing data and learned ensemble similarity network
#'
#' @author Hyundoo Jeong <hdj@inu.ac.kr>
#'
#' @references
#' Accurate single-cell clustering through effective noise reduction over ensemble similarity network
#'

#' @export
siclen_denoise <- function(data, nens = 20, knum = 30, alpha = 0.7, npc = 10, pFeat = 0.05){

  log_scdata <- log2(1+cpm(data))

  myvar_genes <- IdentifyFeatures(log_scdata, pFeat = pFeat)
  var_genes <- myvar_genes$VarGenes
  rname <- rownames(log_scdata)


  ens_net <- make_ensnet(log_scdata, var_genes, K = knum, nPCs = npc, nens=nens)

  knn_ensnet <- drop0(ens_net, tol = 0.5*(nens-1))

  diag(knn_ensnet) <- 0
  idx <- which(colSums(knn_ensnet) == 0)
  knn_ensnet[idx, ] <- ens_net[idx,]
  knn_ensnet[,idx ] <- t(ens_net[idx,])
  rm(ens_net)

  scale_knn <- knn_ensnet %*% Diagonal(x=1/colSums(knn_ensnet))
  scale_knn <- scale_knn %*% scale_knn

  eye <- Diagonal(dim(log_scdata)[2])
  R <- (eye - alpha * scale_knn)
  e <- t((1-alpha) * log_scdata)
  ssp <- solve(R,e)

  ccnt <- 2^ssp - 1
  ccnt <- cpm(t(ccnt))

  results <- list()
  results$corrected_cpm <- ccnt
  results$eknnnet <- knn_ensnet

  return(results)
}







#' SICLEN: Single-cell Clustering based on effective noise reduction through ensenble similarity network
#' @description
#' This function reduces the noise in the single-cell RNA sequencing through a random walk restart over the ensenble similarity network  the number of clusters in the single-cell RNA sequencing data
#'
#' @param path file name to be loaded
#'
#' @return
#' load R objects
#'
#' @author Hyundoo Jeong <hdj@inu.ac.kr>
#'
#' @references
#' Accurate single-cell clustering through effective noise reduction over ensemble similarity network
#'
#' @export
load_obj <- function(path)
{
  env <- new.env()
  nm <- load(path, env)[1]
  env[[nm]]
}






#' SICLEN: Single-cell Clustering based on effective noise reduction through ensenble similarity network
#' @description
#' This function identifies the potential feature genes
#'
#' @param inData single-cell RNA seuqncing data
#' @param  pFeat Percentage of features to be selected
#' @return
#' List of potential feature genes
#'
#' @author Hyundoo Jeong <hdj@inu.ac.kr>
#'
#' @references
#' Accurate single-cell clustering through effective noise reduction over ensemble similarity network
#'
#' @export
IdentifyFeatures <- function(inData = data, pFeat = 0.1){


  rm <- rowMeans(inData)
  rv <- rowVars(inData)

  hind <- which(rm >= median(rm))
  my_rv <- rv[hind]
  names(my_rv) <- names(hind)
  my_vargenes <- names(head(x = sort(my_rv, decreasing = T), ceiling(pFeat*dim(inData)[1])) )

  topK <- names(head(x = sort(my_rv, decreasing = T), ceiling(0.05*dim(inData)[1])) )


  features <- list("VarGenes" = my_vargenes,
                   "topK" = topK

  )
  return(features)
}



#' SICLEN: Single-cell Clustering based on effective noise reduction through ensenble similarity network
#' @description
#' This function constructs the ensemble similarity network based on different cell-to-cell similarity estimations
#'
#' @param dat_mat (M by N) dimensional single-cell RNA seuqncing data
#' @param var_genes The set of potential feature genes
#' @param K The number of nearest neighbors for KNN networks
#' @param nPCs The number of principal components
#' @param nens The number of similarity estimates to construct the ensemble similarity network
#'
#' @return
#' Ensemble similarity network based on the different cell-to-cell similarity estimations
#'
#' @author Hyundoo Jeong <hdj@inu.ac.kr>
#'
#' @references
#' Accurate single-cell clustering through effective noise reduction over ensemble similarity network
#'
#' @export
make_ensnet <- function (dat_mat, var_genes, K = 30, nPCs = 10, nens=20)
{


  Ncells <- dim(dat_mat)[2]
  ngenes <- length(var_genes);

  knn_ens <- Matrix(0, nrow = Ncells, ncol = Ncells, sparse = TRUE)
  sample_rng <- round(runif(n=nens, min = 0.5, max = 0.75), digit = 2)

  for (w in sample_rng){

    sampled_genes <- sample(var_genes, size = ceiling(ngenes*w), replace = F)

    my_pca_log <- prcomp_irlba(t(dat_mat[sampled_genes,]), n = nPCs)


    PCs <- t(my_pca_log$x)
    cor_dist <- cor(PCs, method ='pearson')


    knn_adj <- make.kNNG(1-cor_dist, k = K, symm = F)

    idx <- which(knn_adj > 0)
    th <- quantile(cor_dist[idx], probs = 0.1)
    cor_dist[cor_dist>=th] <- 1
    cor_dist[cor_dist<th] <- 0


    diag(cor_dist) <- 0
    idx <- which(colSums(cor_dist) == 0)
    cor_dist[idx,] <- knn_adj[idx,]
    cor_dist[,idx] <- t(knn_adj[idx,])

    knn_ens <- knn_ens + cor_dist
  }
  return (knn_ens)

}

