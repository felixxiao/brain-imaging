######################################################################

compute.adjacency = function(edges)
{
  .validate.edges(edges)
  
  V = sort(unique(as.numeric(edges$edge.mat)))
  n = length(V)

  edges$edge.mat[,1] = mapvalues(edges$edge.mat[,1], V, 1:n, warn_missing = F)
  edges$edge.mat[,2] = mapvalues(edges$edge.mat[,2], V, 1:n, warn_missing = F)

  assert_that(max(edges$edge.mat) == n)

  # adjacency matrix
  Matrix::sparseMatrix(i = edges$edge.mat[,1],
                       j = edges$edge.mat[,2],
                       x = edges$energy.vec,
                       dims = c(n, n),
                       symmetric = T)
}

compute.laplacian = function(edges)
{
  A = compute.adjacency(edges)
  
  # degree matrix
  D = Matrix::Diagonal(x = colSums(A))

  # Laplacian matrix
  L = as(D - A, 'dgCMatrix')
  assert_that(isSymmetric(L))

  L
}

# Spectral relaxation of ratio-cut minimization
partition.spectral.multiway = function(edges, k, laplacian = NULL,
                                       output.objective = F)
{
  assert_that(k > 1)

  if (is.null(laplacian))
    L = compute.laplacian(edges)

  cat('Computing eigenvectors\n')
  eig = eigs_sym(L, k + 1, 'SM')

  V = eig$vectors[,(k+1):2]
  V = t(apply(V, 1, function(v) v / sqrt(sum(v * v))))

  ## cosine similarity k-means
  # initialize seed centroids
  U = matrix(rnorm(k * k), nrow = k)
  U = apply(U, 2, function(u) u / sqrt(sum(u * u)))
  # k-means
  ITER = 10
  TOL = 10^-4
  sim = rep(NA, length = ITER + 1)
  sim[1] = sum(apply(V %*% U, 1, max))
  for (t in 1:ITER)
  {
    cat(t, ' ')
    part = apply(V %*% U, 1, which.max)
    for (k in 1:k)
    {
      idx = which(part == k)
      if (length(idx) > 1)
        U[,k] = colSums(V[idx,])
      else if (length(idx) == 1)
        U[,k] = V[idx,]
      else
        U[,k] = rnorm(k)
    }
    U = apply(U, 2, function(u) u / sqrt(sum(u * u)))
    sim[t+1] = sum(apply(V %*% U, 1, max))
    if (abs(sim[t+1] - sim[t]) < TOL)
      break
  }

  if (output.objective)
    return(list(part = as.factor(part), obj = sim[1:(t+1)]))
  as.factor(part)
}

# recursive (equal) bipartitioning using Laplacian spectrum
partition.spectral.recursive = function(edges, max.depth)
{
  assert_that(max.depth > 0)
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
