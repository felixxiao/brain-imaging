# Used to compute energy covariance between 1-dimensional variables more
# efficiently than the energy package. See (Szekely 2013) for details.
# x : numeric, univariate data vector
.distance.matrix = function(x)
{
  a = matrix(x, nrow = length(x), ncol = length(x))
  a = abs(a - t(a))
  a_m = matrix(apply(a, 1, mean), nrow = length(x), ncol = length(x))
  a - a_m - t(a_m) + mean(a)
}

# Return a matrix with row i and column j containing the energy correlation
# between column i of dat and column j of dat2 (dat if dat2 is no specified)
.pairwise.ecor = function(dat, dat2 = NULL, na.replace = 1)
{
  n = nrow(dat)
  m = ncol(dat)

  if (! is.null(dat2))
  {
    stopifnot(nrow(dat2) == n)

    dist2 = matrix(0, nrow = n*n, ncol = ncol(dat2))
    for (j in 1:ncol(dat2))
      dist2[,j] = as.numeric(.distance.matrix(dat2[,j]))
  }
  
  dist = matrix(0, nrow = n*n, ncol = ncol(dat))
  for (j in 1:ncol(dat))
    dist[,j] = as.numeric(.distance.matrix(dat[,j]))

  if (is.null(dat2)) dist2 = dist

  ecov = sqrt(crossprod(dist, dist2)) / n
  evar1 = sqrt(colSums(dist * dist)) / n
  evar2 = sqrt(colSums(dist2 * dist2)) / n
  ecor = sweep(ecov, 1, sqrt(evar1), '/')
  ecor = sweep(ecor, 2, sqrt(evar2), '/')
  
  ecor[is.na(ecor)] = na.replace

  ecor
}

# For each component, return the mean and variance of energy correlation
# statistics between all pairs of voxels within the component.
# For computational speed, optionally restrict pairwise ecor computation to a
# subsample of all voxels in each component.
# Arguments
#   dat       : matrix with row as time and column as voxel
#   parcel    : factor of cluster assignments with length
#               equal to the number of voxels
#   subsample : if numeric, pairwise ecor will be computed from at most
#               this number of voxels for each component
#   replace   : logical, default T to subsample with replacement, F otherwise
criterion.within_pairwise_ecor = function(dat, parcel, subsample = 1000,
                                          replace = T, verbose = T)
{
  assert_that(ncol(dat) == length(parcel))

  n = nrow(dat)
  if (! is.numeric(subsample)) subsample = Inf
  
  within.mean = rep(NA, nlevels(parcel))
  within.var  = rep(NA, nlevels(parcel))
  within.n    = rep(NA, nlevels(parcel))
  names(within.var)  = levels(parcel)
  names(within.mean) = levels(parcel)
  names(within.n)    = levels(parcel) 

  if (verbose) cat('Computing Within-Score for Component:')
  for (p in levels(parcel))
  {
    if (verbose) cat(' ', p)

    if (length(which(parcel == p)) <= subsample)
      voxels = which(parcel == p)
    else
      voxels = sample(which(parcel == p), size = subsample, replace = replace)
    
    ecor = .pairwise.ecor(dat[,voxels])
    
    within.mean[p] = mean(ecor)
    within.var[p]  = var(as.numeric(ecor))
    within.n[p] = length(voxels)
  }
  if (verbose) cat('\n')
  
  list(mean = within.mean, var = within.var, n = within.n,
       total.mean = mean(within.mean, na.rm = T))
}

# Arguments
#   dat       : matrix with row as time and column as voxel
#   parcel    : factor of cluster assignments with length
#               equal to the number of voxels
#   subsample : if numeric, pairwise ecor will be computed from at most
#               this number of voxels for each component
#   replace   : logical, default T to subsample with replacement, F otherwise
criterion.between_pairwise_ecor = function(dat, parcel, subsample = 1000,
                                           replace = T, verbose = T)
{
  assert_that(ncol(dat) == length(parcel))

  n = nrow(dat)
  m = nlevels(parcel)
  if (m < 2)
  {
    warning('only one component')
    return()
  }
  if (! is.numeric(subsample)) subsample = Inf
  
  component_pairs = combn(levels(parcel), 2)
  
  between.mean = matrix(NA, nrow = m, ncol = m)
  rownames(between.mean) = levels(parcel)
  colnames(between.mean) = levels(parcel)
  between.var = between.mean
  between.n   = rep(NA, times = nlevels(parcel))
  names(between.n) = levels(parcel)
  
  if (verbose) cat('Computing Between-Score for Components:')
  for (pair in 1:ncol(component_pairs))
  {
    i = component_pairs[1, pair]
    j = component_pairs[2, pair]
    
    if (verbose) cat(' (', i, ',', j, ')', sep = '')
    
    voxelsI = which(parcel == i)
    voxelsJ = which(parcel == j)
    if (length(voxelsI) > subsample)
      voxelsI = sample(voxelsI, subsample, replace)
    if (length(voxelsJ) > subsample)
      voxelsJ = sample(voxelsJ, subsample, replace)
    
    ecor = .pairwise.ecor(dat[,voxelsI], dat[,voxelsJ])

    between.mean[i,j] = between.mean[j,i] = mean(ecor)
    between.var[i,j]  = between.var[j,i]  = var(as.numeric(ecor))
    between.n[i] = length(voxelsI)
    between.n[j] = length(voxelsJ)
  }
  if (verbose) cat('\n')

  list(mean = between.mean, var = between.var, n = between.n,
       total.mean = mean(between.mean, na.rm = T))
}

criterion.adjacent_pairwise_ecor = function(edges, parcel)
{
  adjacent.mean = rep(NA, nlevels(parcel))
  names(adjacent.mean) = levels(parcel)
  adjacent.var  = adjacent.mean

  # list, for each component, the indices of the edges for which the
  #   1st node is in that component
  parcels1 = split(1:nrow(edges$edge.mat), parcel[edges$edge.mat[,1]])
  # ditto for 2nd node
  parcels2 = split(1:nrow(edges$edge.mat), parcel[edges$edge.mat[,2]])

  for (p in levels(parcel))
  {
    ecor = edges$energy.vec[intersect(parcels1[[p]], parcels2[[p]])]

    adjacent.mean[p] = mean(ecor)
    adjacent.var[p] = var(ecor)
  }

  list(mean = adjacent.mean, var = adjacent.var,
       total.mean = mean(adjacent.mean, na.rm = T))
}

criterion.adjacent_pairwise_ecor.slow = function(edges, parcel)
{
  edge.mat = edges$edge.mat
  weights = edges$energy.vec
  assert_that(nrow(edge.mat) == length(weights))

  edge.mat = mapvalues(edge.mat, 1:length(parcel), parcel,
    warn_missing = F)

  idx = (edge.mat[,1] == edge.mat[,2])
  
  adj.score = sapply(split(weights[idx], edge.mat[idx,1]), mean)
  list(mean = adj.score, total.mean = mean(adj.score, na.rm = T))
}

# may be incorrect
criterion.boundary_pairwise_ecor = function(edges, parcel, verbose = F)
{
  levels(parcel) = 1:nlevels(parcel)
  boundary.mean = matrix(NA, nrow = nlevels(parcel), ncol = nlevels(parcel),
                         dimnames = list(levels(parcel), levels(parcel)))
  boundary.var  = boundary.mean
  boundary.n    = boundary.mean
  boundary.n[]  = 0

  # list, for each component, the indices of the edges for which the
  #   1st node is in that component
  parcels1 = split(1:nrow(edges$edge.mat), parcel[edges$edge.mat[,1]])
  # ditto for 2nd node
  parcels2 = split(1:nrow(edges$edge.mat), parcel[edges$edge.mat[,2]])

  component_pairs = combn(levels(parcel), 2)
  for (pair in 1:ncol(component_pairs))
  {
    comp1 = component_pairs[1, pair]
    comp2 = component_pairs[2, pair]
    if (verbose) cat('(', comp1, ', ', comp2, ')')

    ecor = edges$energy.vec[c(intersect(parcels1[[comp1]], parcels2[[comp2]]),
                              intersect(parcels2[[comp1]], parcels2[[comp2]]))]

    if (length(ecor) > 0)
    {
      boundary.mean[comp1, comp2] = boundary.mean[comp2, comp1] = mean(ecor)
      boundary.var[comp1, comp2]  = boundary.var[comp2, comp1]  = var(ecor)
      boundary.n[comp1, comp2]    = boundary.n[comp2, comp1]    = length(ecor)
    }
  }

  total.mean = sum(boundary.mean * boundary.n, na.rm = T) /
               sum(boundary.n, na.rm = T)

  list(mean = boundary.mean, var = boundary.var, n = boundary.n,
       total.mean = total.mean)
}

criterion.multi_boundary_ecor = function(dat, parcel, detailed = F)
{
  assert_that(ncol(dat) == length(parcel))
  parcel.list = split(1:ncol(dat), parcel)
  m = nrow(dat)

  A = sapply(parcel.list, function(cols) {
        a = array(dat[,cols], dim = c(m, length(cols), m))
        a = apply(a - aperm(a, c(3,2,1)), c(1,3),
              function(x) sqrt(sum(x * x)))
        colmeans.a = colMeans(a)
        a = sweep(a, 1, colmeans.a, '-')
        a = sweep(a, 2, colmeans.a, '-')
        as.vector(a + mean(colmeans.a))
      }, simplify = 'array')

  V = crossprod(A)
  sqrt.diag.V = sqrt(diag(V))
  V = sweep(V, 1, sqrt.diag.V, '/')
  V = sweep(V, 2, sqrt.diag.V, '/')
  V = sqrt(V)
  if (detailed) return(V)
  mean(V[upper.tri(V)])
}

# type can be 'ratio', 'normalized'
#   if left unspecified,  
criterion.cut_weight = function(edges, part, type = '')
{
  L = compute.laplacian(edges)

  X = convert.factor.assignment_mat(part, ratioed = type == 'ratio')

  matrix.trace(as.matrix(t(X) %*% L %*% X))
}

criterion.balance = function(part)
{
  sizes = as.vector(table(part))
  size.ratios = sizes / max(sizes)
  mean(size.ratios)
}

criterion.jaggedness = function(edges, part, detailed = F)
{ 
  .validate.edges(edges)
  part.mat = mapvalues(edges$edge.mat, 1:length(part), part,
    warn_missing = F)
  surface = as.vector(part.mat[part.mat[,1] != part.mat[,2],])
  surface = table(surface)
  volume = table(part)
  ratio = (surface)^(3/2) / volume
  if (detailed) return(ratio)
  mean(ratio)
}

criterion.adj_rand_similarity = function(part1, part2)
{

}