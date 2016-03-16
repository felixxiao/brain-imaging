# convert adj.list into a matrix
convert.adjList2edgeMat = function(adj.list, duplicates = F)
{
  v1 = rep(1:length(adj.list), times = lapply(adj.list, length))
  v2 = unlist(adj.list)

  assert_that(length(v1) == length(v2))

  if (! duplicates)
  {
    not.dup = v1 >= v2
    v1 = v1[not.dup]
    v2 = v2[not.dup]
  }

  cbind(v1, v2)
}

# Arguments
#   part  : factor, partitioning assignments of graph A vertices
#   map   : integer, maps vertices from graph A to graph B
#   N     : numeric(1), number of vertices in graph B
# Return
#   factor(N), with partitioning assignments in A mapped to B vertices
#     and NA values for B vertices with no corresponding A vertices
preprocess.map_partition = function(part, map, N)
{
  assert_that(length(part) == length(map))
  assert_that(length(part) <= N)
  
  part.map = factor(rep(NA, times = N), levels = 1:(nlevels(part) + 1))
  part.map[map] = part
  part[is.na(part)] = nlevels(part) + 1
  part.map
}

# Arguments
#   edge.mat  : matrix, each row contains the indices of two vertices
#               joined in an edge
#   map       : numeric, vertex index --> new index (positive integer)
#   inverse   : logical(1), if TRUE, map new index back to original
#   na.rm     : logical(1), if TRUE, any edge containing a vertex that
#               did not map will be removed. If FALSE, edges will be
#               left with NA vertices. Default TRUE.
# Return
#   edge.mat copy with old vertices mapped to new ones
preprocess.map_vertices = function(edge.mat, map, inverse = F,
                                   na.rm = T)
{
  assert_that(all(map > 0))
  if (length(map) > max(edge.mat))
    map = c(map, rep(NA, times = max(edge.mat) - length(map)))
  
  if (inverse)
  {
    assert_that(! any(duplicated(map)))
    map.inv = rep(NA, times = max(edge.mat, map))
    map.inv[map] = 1:length(map)
    map = map.inv
  }

  edge.mat[,1] = mapvalues(edge.mat[,1], 1:length(map), map,
                           warn_missing = F)
  edge.mat[,2] = mapvalues(edge.mat[,2], 1:length(map), map,
                           warn_missing = F)
  if (na.rm)
  {
    na.edges = is.na(edge.mat[,1]) | is.na(edge.mat[,2])
    return(edge.mat[! na.edges,])
  }
  edge.mat
}

# Parameter
#   data  : matrix of fMRI values with voxel columns
# Return
#   $data : matrix of fMRI values with zero-voxel columns removed
#   $map  : numeric, for each voxel column in the new data matrix,
#           which voxel column it corresponded with in the old
preprocess.remove_zero_vertices = function(data)
{
  n.vertices = ncol(data)
  nonzero_vertices = which(apply(data, 2, function(x) any(x != 0)))

  list(data = data[,nonzero_vertices],
       map  = which(1:n.vertices %in% nonzero_vertices))
}

# compute the distance covariance for the edge weights
compute.edgeWeights = function(data, adj.list, func = dcor,
                               edge.mat = NULL, verbose = T, save = F)
{
  if (is.null(edge.mat))
    edge.mat = convert.adjList2edgeMat(adj.list)

  batch.len = ceiling(nrow(edge.mat)/10)
  vec = numeric(nrow(edge.mat))

  # split edge computation into 10 batches
  for (i in 1:10)
  {
    # form the indices we're going to work over
    idx = ((i-1)*batch.len + 1):(min(i*batch.len, nrow(edge.mat)))

    vec[idx] = sapply(idx, function(x){
      func(data[,edge.mat[x,1]], data[,edge.mat[x,2]])
    })
 
    if(verbose) cat('*')
  }
  edges = list(edge.mat = edge.mat, energy.vec = vec)
  
  if (save) save(edges, file = paste0(PATH_SAVE, "edges_", 
                                      Sys.Date(), ".RData"))
  
  edges
}

partition.addedge.unconstrained = function(edges, component.num)
{
  edge.mat = edges$edge.mat
  weights = edges$energy.vec
  assert_that(nrow(edge.mat) == length(weights))
  assert_that(! any(is.na(edge.mat)) & ! any(is.na(weights)))

  edge.mat = edge.mat[order(weights, decreasing = T),]

  n = max(edge.mat)
  g = graph.empty(n, directed = F)
  n.comp = n
  i = 1
  while (n.comp > component.num)
  {
    edges = edge.mat[i:(i + n.comp - component.num - 1),]
    g = add_edges(g, t(edges))
    i = i + n.comp - component.num
    n.comp = components(g)$no
    assert_that(n.comp >= component.num)
    cat(n.comp, ' components\n')
  }

  as.factor(components(g)$membership)
}

# construct a graph by adding edges in descending order of
#   energy dist. stop at specified number of connected components
# data and adj.list are ignored if edges are given
construct.graph <- function(data, adj.list, edges = NULL,
                            component.num = 20,
                            verbose = TRUE, save = TRUE) {
  if (is.null(edges))
  {
    cat('Edges are null\n')
    assert_that(is.numeric(data) & is.matrix(data))
    assert_that(is.list(adj.list))
    
    edges = compute.edgeWeights(data, adj.list, dcor, save = save)
    if(verbose) cat("Edges done computing.")
  }

  #number of vertices
  n = max(edges$edge.mat)
  zero.voxels = which(apply(data, 2, function(x){sum(abs(x))}) == 0) 
  g.base = .construct.graphBase(zero.voxel, edges$edge.mat, n)

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


.construct.graphBase <- function(zero.voxels, edge.mat, n, verbose = FALSE){
  assert_that(is.numeric(n) & length(n) == 1)
  assert_that(is.matrix(edge.mat) & is.numeric(edge.mat))#initialize the graph
 
  #initialize the graph
  g = graph.empty(n, directed = F)
  
  #determine which columns are all 0
  iter = 1

  #do this recursively: find all edge-pairs such that one edge
  #  is not in zero.voxel and one edge is. Add that edge and
  #  move the zero.voxel into the non-zero.voxel list.
  #  Continue until there are no more zero.voxels
  while(length(zero.voxels)>0) {
    bool.mat = apply(edge.mat, 2, function(x){x %in% zero.voxels})

    #find voxels between nonzero and zero
    sum.mat = apply(bool.mat, 1, sum)
    boundary.edgeidx = which(sum.mat == 1)  

    #now for the annoying part: you want to make sure you don't add
    # the same zero.voxel in twice
    boundary.edge = edge.mat[boundary.edgeidx,]
    idx = which(boundary.edge %in% zero.voxels)
    accounted.voxels = boundary.edge[idx]
    duplicated.zerovoxels = duplicated(accounted.voxels)
    #if(sum(duplicated.zerovoxels) > 0) {
      idx = idx[-which(duplicated.zerovoxels == TRUE)]
    #}
 
    #now remove the duplicated zero.voxels
    idx = idx %% nrow(boundary.edge)
    idx[idx == 0] = nrow(boundary.edge)
    boundary.edge = boundary.edge[idx,]

    #add edges to g
    g = add.edges(g, t(boundary.edge))

    #mark zero edges as nonzero
    zero.voxels = zero.voxels[!zero.voxels %in% 
     unique(as.vector(boundary.edge))]

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

