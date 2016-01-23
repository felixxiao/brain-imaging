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
  edges = compute.edgeWeights(mat, neigh$neighbor.list, dcor, save = FALSE,
   verbose = FALSE)

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

  batch.all = unlist(batch.list)
  expect_true(length(batch.all) ==
   length(unique(batch.all)))
  expect_true(5 == 3)

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

test_that("Graph Base is correct", {
  edges = compute.edgeWeights(mat, neigh$neighbor.list, dcor, save = FALSE)

  n = max(edges$edge.mat)
  g.base = .construct.graphBase(mat, edges$edge.mat, n)

  #make sure correct number of vertices
  expect_true(vcount(g.base) == 6^3)

  #make sure all 0-voxels have an edge (degree > 0)
  expect_true(length(which(apply(mat, 2, sum) == 0)) <= 
   length(which(degree(g.base) > 0)))
 
  dimen = c(6,6,6)
  batch = apply(expand.grid(c(1,6), 1:6, 1:6), 1, .convert.3Dto2Dloc,
   dimen = dimen)
  batch = c(apply(expand.grid(1:6, c(1,6), 1:6), 1, .convert.3Dto2Dloc,
   dimen = dimen), batch)
  batch = c(apply(expand.grid(1:6, 1:6, c(1,5,6)), 1, .convert.3Dto2Dloc,
   dimen = dimen), batch)
  expect_true(all(degree(g.base)[batch] > 0))


  edges.gBase = as_edgelist(g.base)
  # make sure there are no self-loops
  bool = apply(edges.gBase, 1, function(x){x[1] == x[2]})
  expect_true(all(!bool))

  #let's check specific cases
  #make sure all the points in X[1, 2:5, 2:5] (zero'd out) are connected
  # to X[2, 2:5, 2:5] respectively
  zero.loc = as.matrix(expand.grid(1, 2:5, 2:4))
  assert_that(all(apply(zero.loc, 1, function(x){all(X[x[1], x[2], x[3], ] 
   == 0)})))
  nonzero.loc = zero.loc
  nonzero.loc[,1] = nonzero.loc[,1] + 1
  assert_that(all(apply(nonzero.loc, 1, function(x){sum(abs(
   X[x[1], x[2], x[3], ])) != 0})))

  check.adjacency <- function(mat1, mat2){
    for(i in 1:nrow(mat1)) {
      idx1 = .convert.3Dto2Dloc(mat1[i,], dimen)
      idx2 = .convert.3Dto2Dloc(mat2[i,], dimen)
  
      #first, make sure zero.idx and nonzero.idx are neighbors
      expect_true(idx2 %in% neigh$neighbor.list[[idx1]])
    
      bool = apply(edges.gBase, 2, function(x){x %in% c(idx1, idx2)})
      bool = apply(bool, 1, sum)
      expect_true(length(which(bool == 2)) >= 1) #can be larger due to duplicates
    }
  }
  check.adjacency(zero.loc, nonzero.loc)

  #similarly, make sure all the points in X[2:5, 2:5, 6] are
  # connected to X[2:5, 2:5, 5], which are in turn connected to
  # X[2:5, 2:5, 4]
  zero.loc = as.matrix(expand.grid(2:5, 2:5, 5))
  nonzero.loc = zero.loc
  nonzero.loc[,3] = nonzero.loc[,3] - 1
  check.adjacency(zero.loc, nonzero.loc)

  zero2.loc = zero.loc
  zero2.loc[,3] = zero2.loc[,3] + 1
  check.adjacency(zero.loc, zero2.loc)

  #check that the points in X[3:4, 3:4, 6] have only one edge coming
  # out, so they should appear in edges.gBase only once
  zero.idx = apply(expand.grid(3:4, 3:4, 6), 1, .convert.3Dto2Dloc, dimen)
  for(i in 1:length(zero.idx)) {
    bool = apply(edges.gBase, 1, function(x){sum(x %in% zero.idx[i])})
    expect_true(sum(bool) <= 2 & sum(bool) >= 1) #handle duplicates for now
  }

  #check the shortest paths matrix
  dmat = shortest.paths(g.base)

  #first check that the edges in the middle (X[3:4, 3:4, 3]) have no
  # edges at all
  loc = apply(expand.grid(3:4, 3:4, 3), 1, .convert.3Dto2Dloc, dimen)
  expect_true(length(unique(as.vector(dmat[loc,]))) == 2) #elements 0 and Inf

  #check that X[3:4, 3:4, 6] are only connected to 2 points (in [,,5] and
  #  [,,4] respectively
  loc = apply(expand.grid(3:4, 3:4, 6), 1, .convert.3Dto2Dloc, dimen)
  for(i in loc) {
    expect_true(length(which(dmat[i,] != Inf)) == 3) 
  }

  #check that the point X[1,1,1] is connected to the point X[2,2,2]
  idx1 = .convert.3Dto2Dloc(c(1,1,1), dimen)
  idx2 = .convert.3Dto2Dloc(c(2,2,2), dimen)
  expect_true(dmat[idx1, idx2] == 3)
})

test_that("Test that construct.graph is correct", {
  g = construct.graph(mat, neigh$neighbor.list, save = FALSE, 
   verbose = FALSE, component.num = 4)

  #make sure all components have more than 2*2*4 members (since
  # by design, for example, all points in X[2:3, 2:3, 2:4] should
  # be in the same group
  comp = components(g)
  tab = table(comp$membership)
  expect_true(all(tab > 2*2*4))

  #check that X[2:3, 2:3, 2:4] are all in the same group
  dimen = c(6,6,6)
  loc = apply(expand.grid(2:3, 2:3, 2:4), 1, .convert.3Dto2Dloc, dimen)
  expect_true(unique(comp$membership[loc]) == 1)
})
