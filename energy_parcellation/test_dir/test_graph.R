#make a dummy dataset that we'll consider
#construction: it's a 6x6x6x10 dataset. The last
# 10 is for time.
#Take a look at the 5x5x6 box: (called X)
# the elements X[,,1], X[,,5], X[,,6],
#              X[1,,], X[6,,],
#              X[,1,], X[,6,] are all 0 (empty)
# the elements X[2:3, 2:3, 2:4] are drawn from one distribution
# the elements X[4:5, 4:5, 2:4] are drawn from another
# the elements X[2:3, 4:5, 2:4] another, and X[4:5, 2:3, 2:4]

set.seed(10)
source("../../source_header.R")
X = array(0, dim = c(6, 6, 6, 10))

#fill in the 4-dimension array X with vec so that 
#  every row along the 4th dimension are correlated
fill.in <- function(dim1, dim2, dim3, vec, noise = 0.05){
  X[dim1, dim2, dim3, ] <<- aperm(apply(X[dim1, dim2, dim3, ,drop=FALSE], 
   c(1:3), function(x){ vec + noise*rnorm(10)}), c(2,3,4,1))
}
#make sure the matrix is adding the way I think it is
vec = rnorm(10)
fill.in(2:3, 2:3, 2:4, vec, 0)
assert_that(all(X[2,2,2,] == X[3,3,4,]))
assert_that(all(X[2,2,2,] != 0))

#now fill in the actual matrix
vec = rnorm(10); fill.in(2:3, 2:3, 2:4, vec)
vec = rnorm(10); fill.in(4:5, 4:5, 2:4, vec)
vec = rnorm(10); fill.in(2:3, 4:5, 2:4, vec)
vec = rnorm(10); fill.in(4:5, 2:3, 2:4, vec)

#now construct the adj.list and 2d matrix
mat = extract.data(X, 1:prod(dim(X)[1:3]))
test_that("extract.data works properly", {
  expect_true(length(which(apply(mat, 2, sum) != 0)) == 4*4*3)
})

#pretend that the entire 6^3 voxels belong in the mask
neigh = extract.neighbors(array(1, dim = c(6,6,6)), 
 pattern = .cross.enumerate())


#now let's start testing
test_that("Make sure that compute.edgeWeights is correct", {
  edges = compute.edgeWeights(mat, neigh$neighbor.list, dcor, save = FALSE)

  #make sure that the dcor in each block is more than the dcor across blocks
  batch.list = list(4)
  batch.list[[1]] = apply(expand.grid(2:3, 2:3, 2:4), 1, .convert.3Dto2Dloc,
   dimen = c(6,6,6))
  batch.list[[2]] = apply(expand.grid(4:5, 4:5, 2:4), 1, .convert.3Dto2Dloc,
   dimen = c(6,6,6))
  batch.list[[3]] = apply(expand.grid(2:3, 4:5, 2:4), 1, .convert.3Dto2Dloc,
   dimen = c(6,6,6))
  batch.list[[4]] = apply(expand.grid(4:5, 2:3, 2:4), 1, .convert.3Dto2Dloc,
   dimen = c(6,6,6))

  expect_true(length(c(batch1, batch2, batch3, batch4)) ==
   length(unique(c(batch1, batch2, batch3, batch4))))

  #find all the energy statistics within a batch and output the median
  output.withinBatch <- function(batch, edges){
    idx1 = edges$edge.mat[,1] %in% batch
    idx2 = edges$edge.mat[,2] %in% batch
    idx = which(apply(cbind(idx1, idx2), 1, sum) == 2)

    median(edges$energy.vec[idx])
  }

  stat.withinBatch = unlist(lapply(batch.list, output.withinBatch, 
   edges = edges))

  #find all the energy statistics across different batches and output the
  #  median
  output.acrossBatch <- function(vec.idx, batch, edges){
    idx1 = edges$edge.mat[,1] %in% batch[[vec.idx[1]]]
    idx2 = edges$edge.mat[,2] %in% batch[[vec.idx[2]]]
    idx = which(apply(cbind(idx1, idx2), 1, sum) == 2)

    idx1 = edges$edge.mat[,1] %in% batch[[vec.idx[2]]]
    idx2 = edges$edge.mat[,2] %in% batch[[vec.idx[1]]]
    idx = c(idx, which(apply(cbind(idx1, idx2), 1, sum) == 2))

    #not every pair of batches have adjacent edges
    if(length(idx)) median(edges$energy.vec[idx]) else 0
  }

  stat.acrossBatch = apply(combn(4,2), 2, output.acrossBatch,
   batch = batch.list, edges = edges)

  expect_true(max(stat.acrossBatch) < min(stat.withinBatch))
})


