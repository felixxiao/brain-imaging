# initialize an undirected graph with no edges and
# size : number of vertices
unionfind = function(size)
{
  list(id = 1:size, num = rep(1, times = size))
}

.unionfind.root = function(uf, i)
{
  while (i != uf$id[i])
  {
    uf$id[i] = uf$id[uf$id[i]]
    i = uf$id[i]
  }
  i
}

# add an edge between vertices i and j
unionfind.union = function(uf, i, j)
{
  rootI = .unionfind.root(uf, i)
  rootJ = .unionfind.root(uf, j)
  if (rootI == rootJ) return(uf)
  if (uf$num[i] > uf$num[j]) 
  {
    uf$id[rootJ] = rootI
    uf$num[rootI] = uf$num[rootI] + uf$num[rootJ]
  }
  else
  {
    uf$id[rootI] = rootJ
    uf$num[rootJ] = uf$num[rootJ] + uf$num[rootI]
  }
  uf
}

# is there a path between vertices i and j?
unionfind.connected = function(uf, i, j)
{
  .unionfind.root(uf, i) == .unionfind.root(uf, j)
}

# how many vertices is i connected to, including itself?
unionfind.component_size = function(uf, i)
{
  uf$num[.unionfind.root(uf, i)]
}
