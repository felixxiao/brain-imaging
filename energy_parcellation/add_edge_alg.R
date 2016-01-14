setwd('C:/Users/Felix/Dropbox/Felix_Kevin_Han-seniorthesis2015-16/data')

#load('edges.RData')
load('ABIDE_50002_matrix_2015-12-07.RData')
load('template_2015-12-07.RData')

library(igraph)
library(energy)


compute.edge.list = function(data, adj.list, func)
{
  v1 = integer()
  v2 = integer()
  w  = numeric()
  for (i in 1:length(adj.list))
  {
    adj = adj.list[[i]]
    for (j in adj)
    {
      w  = c(w, func(data[,i], data[,j]))
      v1 = c(v1, i)
      v2 = c(v2, j)
    }
  }
  return(data.frame(v1 = v1, v2 = v2, w = w))
}

edges = compute.edge.list(dat$dat, template$neighbor.list, dcor)
#save(edges, file = 'edges_list.RData')

load('edges_list.RData')

g = graph.empty(directed = F)
g = add.vertices(g, ncol(dat$dat))
edges = edges[order(edges$w, decreasing = T),]

for (k in 1:nrow(edges))
{
  g = add.edges(g, as.numeric(edges[k,1:2]))

  if (components(g)$no == 20)
    break
}


