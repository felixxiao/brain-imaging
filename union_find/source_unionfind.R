# Usage:
#  uf = UnionFind$new(10)  initialize
#  uf$connected(2,5)       FALSE
#  uf$union(5,2)
#  uf$connected(2,5)       TRUE
#  uf$component_size(5)    2
UnionFind = R6Class('UnionFind',
  public = list(
    id = NA,
    num = NA,
    
    initialize = function(size) {
      self$id = 1:size
      self$num = rep(1, times = size)
    },
    
    # find the component identifier of vertex i
    root = function(i)
    {
      while (i != self$id[i]) {
        self$id[i] = self$id[self$id[i]]
        i = self$id[i]
      }
      i
    },
    
    # add an edge between vertices i and j
    union = function(i, j)
    {
      rootI = self$root(i)
      rootJ = self$root(j)
      if (rootI == rootJ) return()
      if (self$num[rootI] > self$num[rootJ]) {
        self$id[rootJ] = rootI
        self$num[rootI] = self$num[rootI] + self$num[rootJ]
      }
      else {
        self$id[rootI] = rootJ
        self$num[rootJ] = self$num[rootJ] + self$num[rootI]
      }
    },
    
    # is there a path between vertices i and j?
    connected = function(i, j)
    {
      self$root(i) == self$root(j)
    },
    
    # how many vertices is i connected to, including itself?
    component_size = function(i)
    {
      self$num[self$root(i)]
    },
    
    flatten = function()
    {
      roots = sapply(1:length(self$id), self$root)
      self$num = self$num[roots]
      self$id = roots
    }
  )
)

# Find a partition using the Add-Edge algorithm : sequentially add each
# edge between two different components if one of the two below is true
#   1) either component has size < min.size
#   2) the union of both components has size <= max.size
# Arguments
#   edges : list containing [[1]] data.frame listing edges and [[2]] vector
#           of energy correlation weights of those edges
partition.addedge.uf = function(edges,
                                max.size = 15000,
                                min.size = 1000,
                                save.path,
                                verbose = T)
{
  edges$edge.mat = edges$edge.mat[order(edges$energy.vec, decreasing = T),]
  edges$energy.vec = edges$energy.vec[order(edges$energy.vec, decreasing = T)]
  
  n = max(edges$edge.mat)
  uf = UnionFind$new(n)
  
  breaks = seq(0, nrow(edges$edge.mat), by = 25000)
  breaks = c(breaks, nrow(edges$edge.mat))
  
  for (edge.set in 2:length(breaks))
  {
    if (verbose) cat(paste(edge.set - 1, '/', length(breaks) - 1, '\n'))
    for (e in (breaks[edge.set - 1] + 1):breaks[edge.set])
    {
      i = edges$edge.mat[e, 1]
      j = edges$edge.mat[e, 2]
      sizeI = uf$component_size(i)
      sizeJ = uf$component_size(j)
      if ((sizeI < min.size || sizeJ < min.size) ||
          sizeI + sizeJ <= max.size)
        uf$union(i, j)
    }
  }
  
  uf$flatten()
  
  if (! missing(save.path))
    save(uf, file = paste0(PATH_DATA,
                           'edge_connect_',
                           MIN_COMP_SIZE, '_',
                           MAX_COMP_SIZE, '_',
                           DATE,
                           '.RData'))
  
  parcels = as.factor(uf$id)
  levels(parcels) = 1:nlevels(parcels)
  
  parcels
}
