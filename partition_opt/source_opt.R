######################################################################

compute.laplacian = function(edges)
{
  V = sort(unique(as.numeric(edges$edge.mat)))
  n = length(V)

  edges$edge.mat[,1] = mapvalues(edges$edge.mat[,1], V, 1:n, warn_missing = F)
  edges$edge.mat[,2] = mapvalues(edges$edge.mat[,2], V, 1:n, warn_missing = F)

  assert_that(max(edges$edge.mat) == n)

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

# recursive bipartitioning using Laplacian spectrum
partition.spectral = function(edges, max.depth = 1)
{
  as.factor(.partition.spectral(edges$edge.mat, edges$energy.vec, 1, 1, max.depth))
}

.partition.spectral = function(edge.mat, weights, label, depth, max.depth)
{
  V = sort(unique(as.numeric(edge.mat)))
  n = length(V)
  if (depth > max.depth)
    return(rep(label, times = n))

  cat(label, '\n')
  L = compute.laplacian(list(edge.mat = edge.mat, energy.vec = weights))
  fiedler = eigs_sym(L, 2, 'SM')$vectors[,1]

  order = V[order(fiedler)]
  part1 = sort(order[1:floor(n / 2)])
  part2 = sort(order[(floor(n / 2) + 1):n])

  idx1 = edge.mat[,1] %in% part1 & edge.mat[,2] %in% part1
  idx2 = edge.mat[,1] %in% part2 & edge.mat[,2] %in% part2

  partition = integer(n)
  names(partition) = V
  part1 = as.character(part1)
  part2 = as.character(part2)
  partition[part1] = .partition.spectral(edge.mat[idx1,], weights[idx1],
    label, depth + 1, max.depth)
  partition[part2] = .partition.spectral(edge.mat[idx2,], weights[idx2],
    label + 2^(max.depth - depth), depth + 1, max.depth)

  partition
}


write.ampl_data = function(laplacian, num.components, file)
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

compute.factor_qp = function(Z, X_init, lambda, ITER = 10)
{
#  Z = as.matrix(read.csv(csv.file, header = F))
  assert_that(isSymmetric(Z))
  assert_that(class(lambda) %in% c('function', 'numeric'))
  n = nrow(Z)
  X = X_init
  k = ncol(X)
  cost = rep(NA, times = ITER)

  if (class(lambda) == 'function')
    lambda = sapply(1:ITER, lambda)
  else if (length(lambda) == 1)
    lambda = rep(lambda, times = ITER)
  else
    assert_that(length(lambda) == ITER)

  for (t in 1:ITER)
  {
    cat(t, ' ')
    # \| Z - X X^T \|_F^2
    cost[t] = norm(Z - X %*% t(X), 'F')

    # min. 1/2 x^T Q x - p^T x
    p = as.vector(crossprod(X, Z + lambda[t] * diag(n)))
    Q = kronecker(diag(n), crossprod(X)) + lambda[t] * diag(n * k)

    # s.t. A^T x >= b
    A = - kronecker(diag(n), matrix(1, nrow = k))
    A = cbind(A, diag(n*k))
    b = c(rep(-1, times = n), rep(0, times = n*k))

    opt = solve.QP(Q, p, A, b, n)
    X = t(matrix(opt$solution, nrow = k))
  }
  X
}