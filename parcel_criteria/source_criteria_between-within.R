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
    stopifnot(nrow(dat) != n)

    dist2 = matrix(0, nrow = n*n, ncol = ncol(dat2))
    for (j in 1:ncol(dat2))
      dist2[,j] = as.numeric(.distance.matrix(dat2[,j]))
    dist2 = dist2 / diag()
  }
  
  dist = matrix(0, nrow = n*n, ncol = ncol(dat))
  for (j in 1:ncol(dat))
    dist[,j] = as.numeric(.distance.matrix(dat[,j]))

  if (is.null(dat2)) dist2 = dist

  ecov = sqrt(t(dist) %*% dist2)
  evar1 = sqrt(matrix(diag(t(dist) %*% dist),
                      nrow = ncol(dist), ncol = ncol(dist2)))
  evar2 = sqrt(t(matrix(diag(t(dist2) %*% dist2),
                        nrow = ncol(dist2), nrow = ncol(dist))))
  ecor = ecov / sqrt(evar1 * evar2)
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
  n = nrow(dat)
  
  within.mean = rep(NA, nlevels(parcel))
  within.var  = rep(NA, nlevels(parcel))
  names(within.var) = levels(parcel)
  names(within.mean) = levels(parcel)

  if (verbose) cat('Computing Within-Score for Component:')
  for (p in levels(parcel))
  {
    if (verbose) cat(' ', p)

    if (! is.numeric(subsample) | length(which(parcel == p)) <= subsample)
      voxels = which(parcel == p)
    else
      voxels = sample(which(parcel == p), size = subsample, replace = replace)
    
    ecor = .pairwise.ecor(dat[,voxels])
    
    within.mean[p] = mean(ecor)
    within.var[p]  = var(as.numeric(ecor))
  }
  if (verbose) cat('\n')
  
  list(mean = within.mean, var = within.var)
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
  n = nrow(dat)

  component_pairs = combn(levels(parcel), 2)
  m = ncol(component_pairs)

  between.mean = matrix(NA, nrow = m, ncol = m)
  rownames(between.mean) = levels(parcel)
  colnames(between.var)  = levels(parcel)
  between.var  = between.mean
  for (pair in 1:m)
  {
    i = component_pairs[1, pair]
    j = component_pairs[2, pair]
    
    ecor = .pairwise.ecor(dat[,parcel == i], dat[,parcel == j])

    between.mean[i,j] = between.mean[j,i] = mean(ecor)
    between.var[i,j]  = between.var[j,i]  = var(as.numeric(ecor))
  }

  list(mean = between.mean, var = between.var)
}

