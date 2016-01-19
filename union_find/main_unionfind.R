# ------------------ Notes for Felix ------------------------
# UnionFind data structure seems to be working correctly
# Unsure why singleton components still exist
# Need to clean up this script and put into a function

PATH_DATA = 'C:/Users/Felix/DropBox/Felix_Kevin_Han-seniorthesis2015-16/data/'
DATE = Sys.Date()
setwd('~/GitHub/union_find')
source('source_unionfind.R')

MAX_COMP_SIZE = 10000

load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

edges$edge.mat = edges$edge.mat[order(edges$energy.vec, decreasing = T),]
edges$energy.vec = edges$energy.vec[order(edges$energy.vec, decreasing = T)]

n = max(edges$edge.mat)
#uf = unionfind(n)
uf = UnionFind$new(n)

edge.added = logical(nrow(edges$edge.mat))

breaks = seq(0, nrow(edges$edge.mat), by = 25000)
breaks = c(breaks, nrow(edges$edge.mat))
for (edge.set in 2:length(breaks))
{
  cat(paste(edge.set - 1, '/', length(breaks) - 1, '\n'))
  for (e in (breaks[edge.set - 1] + 1):breaks[edge.set])
  {
    i = edges$edge.mat[e, 1]
    j = edges$edge.mat[e, 2]
    if (uf$component_size(i) + uf$component_size(j) <= MAX_COMP_SIZE)
    {
      uf$union(i, j)
      edge.added[e] = T
    }
  }
}

save(uf, edge.added, file = paste0(PATH_DATA,
                                   'edge_connect_',
                                   MAX_COMP_SIZE, '_',
                                   DATE,
                                   '.RData'))

sum(edge.added) / nrow(edges$edge.mat)

g = graph.empty(n, directed = F)
g = add.edges(g, t(edges$edge.mat[edge.added,]))
components = components(g)
sort(table(components$membership), decreasing = T)

singletons = which(components$csize == 1)

.unionfind.root(uf, singletons[2])

i = 1
edges$edge.mat[edges$edge.mat[,1] == singletons[i],]
e = which(edges$edge.mat[,1] == singletons[i])
edge.added[e]
edges$edge.mat[e,]
