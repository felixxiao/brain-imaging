# Implementation of (Kim 2011) Fast Nonnegative Matrix Factorization
# block pivoting method, specialized for SymNMF
# Arguments:
#   A     : matrix, symmetric, possibly sparse (n x n)
#   W     : matrix (n x k)
#   alpha : positive scalar
compute.nnls = function(A, W, alpha)
{
  n = nrow(A)
  k = ncol(W)
  assert_that(Matrix::isSymmetric(A))
  assert_that(nrow(W) == n)

# just for easy comparison with notation in (Kim 2011)
#  p = n + k
#  q = k
#  r = n

  CTC = crossprod(W) + alpha * diag(k)
  CTB = Matrix::crossprod(W, A) + alpha * t(W)

  F = matrix(FALSE, k, n)
  G = ! F

  X = matrix(0, k, n)
  Y = - CTB

  compute.XY = function(F)
  {
    # find all j in {1 ... n} that have identical F[,j]
    col_group = split(1:n, crossprod(F, 2^(0:(k-1))))
    # col_group is a list, for each column group, the j indices in it

    X = matrix(0, k, n)
    Y = X
    for (j in col_group)
    {
      f = F[,min(j)]
      g = !f
      if (sum(f) > 0) X[f,j] = as.matrix(solve(CTC[f,f]) %*% CTB[f,j])
      Y[g,j] = as.matrix(CTC[g,f, drop = FALSE] %*% X[f,j, drop = FALSE] - CTB[g,j])
    }
    list(X = X, Y = Y)
  }

  while (any(X < 0) | any(Y < 0))
  {
    a = integer(n) + 3
    b = integer(n) + k + 1

    V = (X < 0) | (Y < 0)

    j = colSums(V) < b
    b[j] = colSums(V)[j]
    a[j] = 3

    j = (colSums(V) >= b) & (a >= 1)
    a[j] = a[j] - 1

    j = (colSums(V) >= b) & (a == 0)
    if (sum(j) > 0) cat(sum(j), ' ')
    V[,j] = apply(V[,j], 2, function(v_j) {
        v = logical(k)
        v[which.max(v_j)] = TRUE
        v
      })

    assert_that(all(F[as.matrix(V)] | G[as.matrix(V)]))

    old.F = F
    F = matrix(FALSE, k, n)
    F[as.matrix(old.F & !V) | as.matrix(V & G)] = TRUE
    G = !F

    XY = compute.XY(F)
    X = XY$X
    Y = XY$Y
  }     

  X
}

# Implementation of (Kuang 2015) SymNMF algorithm for graph partitioning
# Arguments
#   A : 
partition.sym_nmf.adj = function(A, k, alpha, growth = .01,
                                 tol = 10^(-2), iter = 100)
{
  assert_that(Matrix::isSymmetric(A))
  assert_that(alpha > 0) 
  n = nrow(A)
  
  H = matrix(runif(n * k), n, k)

  obj = function(W, H)
  {
    X = (W + H) / 2
    norm(A - X %*% t(X), 'f')^2
  }

  g = c()

  for (t in 1:iter)
  {
    cat(t, ' ')
    W = H
    H = t(compute.nnls(A, W, alpha))
    g = c(g, obj(W, H))
    alpha = alpha * (1 + growth)
  }

  list(W = W, H = H, obj = g)
}

partition.sym_nmf = function(edges, num_components)
{
  A = compute.adjacency(edges)
  k = num_components
  X = partition.sym_nmf.adj(A, k, 1)$H
  as.factor(apply(X, 1, which.max))
}

partition.sym_nmf.python = function(edges, num_components)
{
  .validate.edges(edges)
  setwd('partition_nmf')
  write.table(cbind(edges$edge.mat, edges$energy.vec),
              file = 'edges.txt', row.names = F, col.names = F)
  system(paste('source_symnmf.py edges.txt', num_components))
  X = read.csv('X_symnmf.csv', header = F)
  setwd('..')
  as.factor(apply(X, 1, which.max))
}