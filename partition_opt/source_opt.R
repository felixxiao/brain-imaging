######################################################################

compute.laplacian = function(edges)
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

  L
}

# only do bipartition for now
partition.spectral = function(edges, balanced = T)
{
  L = compute.laplacian(edges)

  ret = eigs_sym(L, 2, 'SM')

  fiedler = ret$vectors[,2]

  partition = factor(fiedler <= quantile(fiedler, 0.5), labels = c(1,2))
  
  partition
}

write.ampl_data = function(laplacian, num.components, file, verbose = T)
{
  L = laplacian
  k = num.components
  n = nrow(L)

  sink(file)

  cat('param n := ', n, ';\n\n', sep = '')
  cat('param k := ', k, ';\n\n', sep = '')

  cat('set V := ')
  cat(1:n)
  cat(';\n\n')

  cat('param L : ')
  cat(sprintf('%8d', 1:n))
  cat(' :=\n')
  for (i in 1:n)
  {
    cat(sprintf('%9d ', i))
    cat(sprintf('%.6f', L[i,]))
    cat('\n')
  }
  cat(';\n\n')

  sink()
}

# needs work
compute.factor_01.julia = function(csv.file, num.components, iter = 5)
{
  julia.file = 'partition_opt/source_factor.jl'
  shell(paste('julia', julia.file, csv.file, num.components, iter))

  X = read.table('partition_opt/X.csv', sep = ',')
  parcel = rep(NA, times = nrow(X))
  for (j in 1:ncol(X)) parcel[which(X[,j])] = j

  as.factor(parcel)
}
