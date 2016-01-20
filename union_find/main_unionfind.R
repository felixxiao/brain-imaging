# ------------------ Notes for Felix ------------------------
# UnionFind data structure seems to be working correctly
# Unsure why singleton components still exist
# Need to clean up this script and put into a function

PATH_DATA = 'C:/Users/Felix/DropBox/Felix_Kevin_Han-seniorthesis2015-16/data/'
DATE = Sys.Date()
setwd('~/GitHub/union_find')
source('source_unionfind.R')

MAX_COMP_SIZE = 15000
MIN_COMP_SIZE =  1000

load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

edges$edge.mat = edges$edge.mat[order(edges$energy.vec, decreasing = T),]
edges$energy.vec = edges$energy.vec[order(edges$energy.vec, decreasing = T)]

n = max(edges$edge.mat)

uf = UnionFind$new(n)

breaks = seq(0, nrow(edges$edge.mat), by = 25000)
breaks = c(breaks, nrow(edges$edge.mat))
for (edge.set in 2:length(breaks))
{
  cat(paste(edge.set - 1, '/', length(breaks) - 1, '\n'))
  for (e in (breaks[edge.set - 1] + 1):breaks[edge.set])
  {
    i = edges$edge.mat[e, 1]
    j = edges$edge.mat[e, 2]
    sizeI = uf$component_size(i)
    sizeJ = uf$component_size(j)
    if ((sizeI < MIN_COMP_SIZE || sizeJ < MIN_COMP_SIZE) ||
        sizeI + sizeJ <= MAX_COMP_SIZE)
      uf$union(i, j)
  }
}

save(uf, edge.added, file = paste0(PATH_DATA,
                                   'edge_connect_',
                                   MIN_COMP_SIZE, '_',
                                   MAX_COMP_SIZE, '_',
                                   DATE,
                                   '.RData'))

uf$flatten()
