#onvert adj.list into a matrix
#TODO: There are duplicates in this list
# Felix: please keep this version with duplicates -- my code uses them
convert.adjList2edgeMat <- function(adj.list) {
  v1 = rep(1:length(adj.list), times = lapply(adj.list, length))
  v2 = unlist(adj.list)

  assert_that(length(v1) == length(v2))

  cbind(v1, v2)
}

#compute the distance covariance for the edge weights
compute.edgeWeights = function(data, adj.list, func = dcor,
                               verbose = FALSE, save = TRUE) {
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
                            verbose = TRUE, save = TRUE) {
  assert_that(is.numeric(data) & is.matrix(data))
  assert_that(is.list(adj.list))

  if(is.null(edges)) {
    edges = compute.edgeWeights(data, adj.list, dcor, save = save)
    if(verbose) cat("Edges done computing.")
  }

  #number of vertices
  n = max(edges$edge.mat)
  g.base = .construct.graphBase(data, edges$edge.mat, n)

  #for each voxel, add the largest associated edge
  #g.base = add.initalEdges(g.base, data, edges) 

  idx = order(edges$energy.vec, decreasing = T)
  edge.mat = edges$edge.mat[idx,]

  #let's do our binary search
  upper = nrow(edge.mat)
  lower = 1
  compute.mid <- function(x,y){floor((x+y)/2)}

  while(TRUE) {
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

#add the largest edge for each voxel
#TO BE FINISHED
.add.initialEdges <- function(g, data, edges){
  idx = order(edges$energy.vec, decreasing = T)
  edge.mat = edges$edge.mat[idx,]

  
}


.construct.graphBase <- function(data, edge.mat, n, verbose = FALSE){
  assert_that(is.numeric(n) & length(n) == 1)
  assert_that(is.matrix(edge.mat) & is.numeric(edge.mat))#initialize the graph
  assert_that(is.matrix(data) & is.numeric(data))

  #initialize the graph
  g = graph.empty(n, directed = F)
  
  #determine which columns are all 0
  zero.voxels = which(apply(data, 2, function(x){sum(abs(x))}) == 0)
  iter = 1

  #do this recursively: find all edge-pairs such that one edge
  #  is not in zero.voxel and one edge is. Add that edge and
  #  move the zero.voxel into the non-zero.voxel list.
  #  Continue until there are no more zero.voxels
  while(length(zero.voxels)>0) {
    bool.mat = apply(edge.mat, 2, function(x){x %in% zero.voxels})

    #find voxels between nonzero and zero
    sum.mat = apply(bool.mat, 1, sum)
    boundary.voxels = which(sum.mat == 1)  

    #add edges to g
    g = add.edges(g, t(edge.mat[boundary.voxels,]))

    #mark zero edges as nonzero
    zero.voxels = zero.voxels[!zero.voxels %in% 
     unique(as.vector(edge.mat[boundary.voxels,]))]

    #remove extra edges
    edge.mat = edge.mat[-which(sum.mat == 0),, drop = FALSE]

    if(verbose) {
      print(iter)
      iter = iter + 1
      print(paste0("Length of zero.voxels: ", length(zero.voxels)))
      print(paste0("Edge count: ", ecount(g)))
      print(paste0("Dim of edges: ", nrow(edge.mat)))
    }
  }

  g
}

