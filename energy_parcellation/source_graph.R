#onvert adj.list into a matrix
convert.adjList2edgeMat <- function(adj.list)
{
  v1 = rep(1:length(adj.list), times = lapply(adj.list, length))
  v2 = unlist(adj.list)

  assert_that(length(v1) == length(v2))

  cbind(v1, v2)
}

#compute the distance covariance for the edge weights
compute.edgeWeights = function(data, adj.list, func = dcor,
                               verbose = FALSE, save = TRUE)
{
  edge.mat = convert.adjList2edgeMat(adj.list)

  batch.len = ceiling(nrow(edge.mat)/10)
  vec = numeric(nrow(edge.mat))

  #split edge computation into 10 batches
  for(i in 1:10){
    #form the indices we're going to work over
    idx = ((i-1)*batch.len+1):(min(i*batch.len, nrow(edge.mat)))

    vec[idx] = sapply(idx, function(x){
      func(data[,edge.mat[x,1]], data[,edge.mat[x,2]])
    })
 
    if(verbose) cat('*')
  }
  edges = list(edge.mat = edge.mat, energy.vec = vec)
  
  if(save) save(edges, file = paste0(PATH_SAVE, "edges_", 
                                     DATE, ".RData"))
  
  edges
}

# construct a graph by adding edges in descending order of
#   energy dist. stop at specified number of connected components
# data and adj.list are ignored if edges are given
construct.graph <- function(data, adj.list, edges = NULL,
                            component.num = 20,
                            verbose = TRUE, save = TRUE)
{
  if(is.null(edges))
  {
    edges = compute.edgeWeights(data, adj.list, dcor, save = save)
    if(verbose) cat("Edges done computing.")
  }
  
  idx = order(edges$energy.vec, decreasing = T)
  edge.mat = edges$edge.mat[idx,]

  # number of vertices
  n = max(edges)
  
  g.base = graph.empty(n, direct = F)

  
  #let's do our binary search
  upper = nrow(edge.mat)
  lower = 1
  compute.mid <- function(x,y){floor((x+y)/2)}

  while(TRUE){
    mid = compute.mid(lower, upper)
    g = add.edges(g.base, t(edge.mat[1:mid,]))

    #add edges until you have component.num components. 
    # too few comps -> add less edges
    # too many comps -> add more edges
    current.comp = components(g)$no
    if(current.comp == component.num) break()
    if(components(g)$no > component.num) lower = mid else upper = mid 

    if(verbose) print(paste0("Current component number: ", current.comp))
  }
  
  g
}

