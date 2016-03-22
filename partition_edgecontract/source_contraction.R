# TODO: add comments

ContractibleGraph = R6Class('ContractibleGraph',
  public = list(
    # list, for each component, the vertices that constitute it
    vertices = NA,
    # list, for each component, list indexed by adjacent components
    # of connecting edge weights as numeric vectors
    edges = NA,
    # number of vertices
    n = NA,
    
    # Assumptions
    # - edges has NO duplicates
    # - in edges$edge.mat, vertices are numbered 1:n
    initialize = function(edges, verbose = F)
    {
      if (is.null(edges)) return()
      edge.mat = edges$edge.mat
      self$n = max(edge.mat)
      energy.vec = edges$energy.vec
      assert_that(length(energy.vec) == nrow(edge.mat))
      
      n = max(edge.mat)
      self$vertices = as.list(1:n)
      names(self$vertices) = as.character(1:n)
      self$edges = lapply(1:n, function(x) list())
      names(self$edges) = as.character(1:n)
      
      breaks = seq(0, nrow(edge.mat), by = 100000)
      breaks = c(breaks, nrow(edge.mat))
      for (i in 1:(length(breaks) - 1))
      {
        if (verbose) cat(i, '/', length(breaks) - 1, ' ')
        for (e in (breaks[i] + 1):breaks[i+1])
        {
          a = as.character(edge.mat[e,1])
          b = as.character(edge.mat[e,2])
          self$edges[[a]][[b]] = self$edges[[b]][[a]] = energy.vec[e]
        }
      }
    },
    
    get_component_names = function()
    {
      names(self$vertices)
    },
    
    get_edges = function(comp)
    {
      c = as.character(comp)
      if (length(self$edges[[c]]) == 0) return(NULL)
      weights = sapply(self$edges[[c]], mean)
      names(weights) = names(self$edges[[c]])
      weights
    },
    
    get_size = function(comp)
    {
      c = as.character(comp)
      length(self$vertices[[c]])
    },
    
    get_adjacency_matrix = function()
    {
    	comp = names(self$vertices)
    	n = length(comp)

      map = 1:n
      names(map) = comp

      A = matrix(0, nrow = n, ncol = n)
      for (c in 1:n)
      {
        edges = self$get_edges(comp[c])
        A[c, map[names(edges)]] = edges
      }
      A
    },

    contract_components = function(comp1, comp2)
    {
      a = as.character(comp1)
      b = as.character(comp2)

      if (! all(c(a, b) %in% self$get_component_names()))
        error(a, ' or ', b, ' is not a component')

      if (is.null(self$edges[[a]][[b]]) | is.null(self$edges[[b]][[a]]))
      {
        assert_that(is.null(self$edges[[b]][[a]]) &
                    is.null(self$edges[[a]][[b]]))
        warning(comp1, ' and ', comp2, ' are not connected')
        return()
      }
      
      assert_that(length(intersect(self$vertices[[a]],
                                   self$vertices[[b]])) == 0)
      self$vertices[[a]] = c(self$vertices[[a]], self$vertices[[b]])
      
      for (b.adj in names(self$edges[[b]]))
      {
        self$edges[[a]][[b.adj]] = c(self$edges[[a]][[b.adj]],
                                     self$edges[[b]][[b.adj]])
        self$edges[[b.adj]][[a]] = self$edges[[a]][[b.adj]]
        self$edges[[b.adj]][[b]] = NULL
      }
      self$edges[[a]][[a]] = NULL

      self$vertices[[b]] = NULL
      self$edges[[b]] = NULL
      
      a
    },
    
    get_vertex_components = function(as.factor = T)
    {
      components = rep(NA, times = self$n)
      for (c in names(self$vertices))
        components[self$vertices[[c]]] = c
      
      if (as.factor)
      {
        components = as.factor(components)
        levels(components) = 1:length(self$vertices)
      }
      
      components 
    }
  )
)

partition.contractedge.python = function(num_components, edges)
{
  write.table(edges$edge.mat, 'partition_edgecontract/edge_mat.csv',
              sep = ',', eol = '\n', row.names = F, col.names = F)
  write.table(edges$energy.vec, 'partition_edgecontract/weights.csv',
              eol = '\n', row.names = F, col.names = F)
  setwd('partition_edgecontract')
  command = 'python source_contraction_singleedge.py'
  system(paste(command, paste(num_components, collapse = ' ')))
  setwd('..')
  
  partitions = lapply(strsplit(
    readLines('partition_edgecontract/partition.csv'), ','),
    as.integer)
  assert_that(all(num_components == partitions[[1]]))
  V = partitions[[2]]
  partitions = lapply(partitions[3:length(partitions)], function(f) {
    part = rep(NA, length = length(map_nz))
    part[V] = f
    as.factor(part)
  })
  list(partitions = partitions, num_components = num_components)
}

# data is either edges or a Contractible Graph
partition.contractedge = function(num.components, data,
                                  return.graph = F, filename = NULL,
                                  verbose = T)
{
  if (class(data)[1] != 'ContractibleGraph')
  {
    assert_that(num.components < max(data$edge.mat))
    if (verbose) cat('Initializing graph ')
    cg = ContractibleGraph$new(data, verbose)
    if (verbose) cat('\n')
  }
  else
    cg = data
  
  if (verbose) cat('Computing initial order\n')
  cg.components = cg$get_component_names()
  n = length(cg.components)
  assert_that(num.components < n)

  # PQ not updated immediately when priorities change (from contraction)
  pq = MinPriorityQueue$new()
  # vectors updated immediately
  priority = c()
  endpoint = c()

  # updates the vectors; optionally updates the PQ
  compute.priority = function(c, insert = F)
  {
    c = as.character(c)
    c.edges = cg$get_edges(c)
    if (length(c.edges) > 0)
    {
      priority[c] <<- cg$get_size(c) - max(c.edges)
      endpoint[c] <<- names(which.max(c.edges))
      if (insert) pq$insert(c, priority[c])
    }
    else # c is either disconnected or no longer exists
      priority[c] <<- Inf
  }
  
  # initialize vectors and PQ
  lapply(cg.components, compute.priority, insert = T)

  if (verbose) cat('Contracting edges; no. components =')
  i = n
  while (i > num.components)
  {
    if (i %% 10000 == 0)
      save(pq, cg, priority, file = paste0('part_ce_', DATE, '.RData'))
    
    pq.min = pq$remove_min()
    a = pq.min$item
    if (is.infinite(priority[a])) {}           # a does not exist
    else if (priority[a] != pq.min$priority)   # a doesn't have this
      pq$insert(a, priority[a])                #   link weight
    else
    {
      b = endpoint[a]
      assert_that(! is.infinite(priority[b]))
      cg$contract_components(a, b)
      i = i - 1
      
      # update priority values of neighbors in the vector, not the PQ
      lapply(names(cg$get_edges(a)), compute.priority)
      # re-insert a back into PQs
      compute.priority(a, insert = T)
      priority[b] = Inf

      if (verbose & i %% 10000 == 0) cat(' ', i)
    }
  }
  if (verbose) cat('\n')
  
  if (! is.null(filename))
    save(cg, file = paste0(PATH_DATA, 'results/cg_',
                           filename, '_',
                           num.components, '_',
                           DATE,
                           '.RData'))

  if (return.graph) return(cg)
  cg$get_vertex_components()
}

"
edges = list()
edges$edge.mat = t(matrix(c(1, 2,
                            1, 3,
                            2, 4,
                            3, 4),
                          nrow = 2))
edges$energy.vec = c(0.7, 0.8, 0.5, 0.9)
#        0.7
#    (1) --- (2)
# 0.8 |       | 0.5
#    (3) --- (4)
#        0.9

cg = ContractibleGraph$new(edges)
cg$get_edges(1)              # 0.7 0.8
cg$get_size(4)               # 1
cg$contract_components(1, 4) # warning, not connected
cg$contract_components(1, 2) # new component 1
cg$get_size(1)               # 2
cg$get_edges(1)              # 0.8 0.5
cg$contract_components(3, 4) # new component 3
cg$get_edges(3)              # (0.8 + 0.5) / 2 = 0.65

cg$contract_components(1, 3)
cg$get_size(1)               # 4

partition.contractedge(edges, num.components = 2) # 1 2 1 1
"
