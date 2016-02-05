library(rARPACK)
library(Matrix)
######################################################################

# only do bipartition for now
partition.spectral = function(edges, balanced = T)
{
  n = max(edges$edge.mat)
  
  # adjacency matrix
  A = Matrix::sparseMatrix(i = edges$edge.mat[,1],
                           j = edges$edge.mat[,2],
                           x = edges$energy.vec,
                           dims = c(n, n),
                           symmetric = T)
  
  # degree matrix
  D = Matrix::Diagonal(x = colSums(A))

  # Laplacian matrix
  L = as(D - A, 'dgCMatrix')
  assert_that(isSymmetric(L))

  ret = eigs_sym(L, 2, 'SM')

  fiedler = ret$vectors[,2]

  partition = factor(fiedler <= quantile(fiedler, 0.5), labels = c(1,2))
  
  partition
}

