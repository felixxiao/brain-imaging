library(R6)

"
.unionfind.RC = setRefClass('UnionFind',
  fields = list(id = 'numeric', num = 'numeric'),
  
  methods = list(
    # find the component identifier of vertex i
    root = function(i)
    {
      while (i != id[i]) {
        id[i] <<- id[id[i]]
        i = id[i]
      }
      i
    },
    
    # add an edge between vertices i and j
    union = function(i, j)
    {
      rootI = root(i)
      rootJ = root(j)
      if (rootI == rootJ) return()
      if (num[rootI] > num[rootJ]) {
        id[rootJ] <<- rootI
        num[rootI] <<- num[rootI] + num[rootJ]
      }
      else {
        id[rootI] <<- rootJ
        num[rootJ] <<- num[rootJ] + num[rootI]
      }
    },
    
    # is there a path between vertices i and j?
    connected = function(i, j)
    {
      root(i) == root(j)
    },
    
    # how many vertices is i connected to, including itself?
    component_size = function(i)
    {
      num[root(i)]
    }
  )
)
"

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

"
# initialize an undirected graph with no edges and
# size : number of vertices
unionfind = function(size)
{
  uf = .unionfind()
  
  uf$id = 1:size
  uf$num = rep(1, times = size)
  
  uf
}
"

# uf = UnionFind$new(10)
# uf$connected(2,5)
# uf$union(5,2)
# uf$connected(2,5)
# uf$component_size(5)
