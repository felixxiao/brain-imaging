source('header_global.R')
library(plyr)
library(gtools)

neighbor.indices = function(dist_sqrd)
{
  n = ceiling(sqrt(dist_sqrd)) + 1
  A = array(0, dim = c(n, n, n))
  for (i in 1:n)
    for (j in 1:n)
      for (k in 1:n)
        A[i,j,k] = (i-1)^2 + (j-1)^2 + (k-1)^2
  idx.nonneg = arrayInd(which(A == dist_sqrd), dim(A)) - 1
  if (nrow(idx.nonneg) == 0)
  {
    warning(dist_sqrd, ' has no valid vertex pairs')
    return()
  }
  signs = gtools::permutations(2, 3, c(1,-1), repeats.allowed = T)
  idx = mapply(kronecker, split(signs, col(signs)),
               split(idx.nonneg, col(idx.nonneg)))
  unname(unique(idx, MARGIN = 1))
}

edgemat = function(dist_sqrd, loc, dimen, mask)
{
  idx = neighbor.indices(dist_sqrd)
  if (is.null(idx)) return()
  edge.mat = lapply(1:ncol(loc), function(j)
    {
      neighbors.3D = sweep(idx, 2, loc[,j], '+')
      neighbors.3D = neighbors.3D[apply(neighbors.3D, 1,
                                        function(row) all(row > 0)),]
      neighbors.1D = apply(neighbors.3D, 1,
                           .convert.3Dto2Dloc, dimen = dimen)
      cbind(mask[j], neighbors.1D, deparse.level = 0)
    })
  edge.mat = do.call(rbind, edge.mat)
  edge.mat = t(apply(edge.mat, 1, sort))
  edge.mat = unique(edge.mat, MARGIN = 1)
  edge.mat = edge.mat[edge.mat[,1] %in% mask & edge.mat[,2] %in% mask,]

  plyr::mapvalues(edge.mat, mask, 1:length(mask), warn_missing = F)
}

