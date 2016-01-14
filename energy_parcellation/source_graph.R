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

construct.graph <- function(data, adj.list, component.num = 20,
 verbose = TRUE){

  edges = compute.edge.list(data, adj.list, dcor)
  save(edges, file = paste0(PATH_DATA, "edges_", DATE, ".RData"))
  if(verbose) cat("Edges done computing.")

  g = graph.empty(direct = F)
  g = add.vertices(g, ncol(dat$dat))
  edges = edges[order(edges$w, decreasing = T),]

  for (k in 1:nrow(edges))
  {
    g = add.edges(g, as.numeric(edges[k,1:2]))

    if (components(g)$no == component.num) break()
  }
}

