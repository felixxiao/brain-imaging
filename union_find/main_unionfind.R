# ------------------ Notes for Felix ------------------------
# UnionFind data structure seems to be working correctly
# Unsure why singleton components still exist
# Need to clean up this script and put into a function


source('source_unionfind.R')

MAX_COMP_SIZE = 10000

load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

partition.addedge.
edges$edge.mat = edges$edge.mat[order(edges$energy.vec, decreasing = T),]
edges$energy.vec = edges$energy.vec[order(edges$energy.vec, decreasing = T)]

n = max(edges$edge.mat)
uf = unionfind(n)

edge.added = logical(nrow(edges$edge.mat))

breaks = seq(0, nrow(edges$edge.mat), by = 25000)
breaks = c(breaks, nrow(edges$edge.mat))
for (edge.set in 2:length(breaks))
{
  cat(paste(edge.set - 1, '/', length(breaks) - 1, '\n'))
  for (e in (breaks[edge.set - 1] + 1):breaks[edge.set]) #nrow(edges$edge.mat))
  {
    i = edges$edge.mat[e, 1]
    j = edges$edge.mat[e, 2]
    if (unionfind.component_size(uf, i) +
        unionfind.component_size(uf, j) <= MAX_COMP_SIZE)
    {
      uf = unionfind.union(uf, i, j)
      edge.added[e] = T
    }
  }
}

#uf = unionfind.flatten(uf)
#roots = sapply(100000:101000, .unionfind.root, uf = uf)

nrow(edges$edge.mat)
sum(edge.added)

g = graph.empty(n, directed = F)
g = add.edges(g, t(edges$edge.mat[edge.added,]))
components = components(g)
sort(table(components$membership), decreasing = T)

tmp = which(components$membership == which.max(table(components$membership)))

singletons = which(components$csize == 1)

i = 1
edges$edge.mat[edges$edge.mat[,1] == singletons[i],]
e = which(edges$edge.mat[,1] == singletons[i] & edges$edge.mat[,2] == 42)
