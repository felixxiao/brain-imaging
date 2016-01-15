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
  for(i in 1:10){
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
 verbose = TRUE, save = TRUE){

  edges = compute.edgeWeights(data, adj.list, dcor)
  if(save) save(edges, file = paste0(PATH_DATA, "edges_", 
   DATE, ".RData"))
  if(verbose) cat("Edges done computing.")

  idx = order(edges$energy.vec, decreasing = TRUE)
  adj.mat = edges$adj.mat[idx,]

  g.base = graph.empty(ncol(dat$dat), direct = F)

  #let's do our binary search
  upper = nrow(adj.mat)
  lower = 1
  compute.mid <- function(x,y){floor((x+y)/2)}

  while(TRUE){
    mid = compute.mid(lower, upper)
    g = add.edges(g.base, t(adj.mat[1:mid,]))

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

