#onvert adj.list into a matrix
convert.edgelist2mat <- function(adj.list){
  vec1 = rep(1:length(adj.list), times = lapply(adj.list, length))
  vec2 = unlist(adj.list)

  assert_that(length(vec1) == length(vec2))

  cbind(vec1, vec2)
}

#compute the distance covariance for the edge weights
compute.edgeWeights = function(data, adj.list, func, verbose = FALSE)
{
  adj.mat = convert.edgelist2mat(adj.list)

  batch.len = ceiling(nrow(adj.mat)/10)
  vec = numeric(nrow(adj.mat))

  #split edge computation into 10 batches
  for(i in 2:10){
    #form the indices we're going to work over
    idx = ((i-1)*batch.len+1):(min(i*batch.len, nrow(adj.mat)))

    vec[idx] = sapply(idx, function(x){
      dcor(data[,adj.mat[x,1]], data[,adj.mat[x,2]])
    })
 
    if(verbose) cat('*')
  }

  list(adj.mat = adj.mat, energy.vec = vec)
}

#construct a graph, edges weighed by energy dist. stop at
#  specified number of connected components
construct.graph <- function(data, adj.list, component.num = 20,
 verbose = TRUE){

  edges = compute.edgeWeights(data, adj.list, dcor)
  save(edges, file = paste0(PATH_DATA, "edges_", DATE, ".RData"))
  if(verbose) cat("Edges done computing.")

  idx = order(edges$energy.vec, decreasing = TRUE)
  adj.mat = edges$adj.mat[idx,]

  g = graph.empty(direct = F)
  g = add.vertices(g, ncol(dat$dat))
 
  for (k in 1:nrow(adj.mat)){
    g = add.edges(g, adj.mat[k,1:2])

    current.comp = components(g)$no 
    if (current.comp == component.num) break()

    if(verbose && k %% floor(nrow(adj.mat)/100) == 0) 
      print(current.comp)
  }
}

