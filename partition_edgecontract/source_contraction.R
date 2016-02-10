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
      n.comp = length(self$vertices)
      A = matrix(0, nrow = n.comp, ncol = n.comp)
      for (c in 1:n.comp)
      {
        edges = self$get_edges(c)
        A[c, as.integer(names(edges))] = edges
      }
      A
    },

    contract_components = function(comp1, comp2)
    {
      a = as.character(comp1)
      b = as.character(comp2)

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

# data is either edges or a Contractible Graph
partition.contractedge = function(num.components, data,
                                  return.graph = F, save.path, verbose = T)
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
  priority = rep(NA, times = n)
  names(priority) = cg.components
  for (c in cg.components)
    priority[c] = 1 - max(cg$get_edges(c))
  
  if (verbose) cat('Contracting edges; no. components =')
  for (i in n:(num.components + 1))
  {
    if (verbose & i %% 10000 == 0) cat(' ', i)
    a = names(which.min(priority))
    b = names(which.max(cg$get_edges(a)))
    cg$contract_components(a, b)

    priority[a] = cg$get_size(a) - max(cg$get_edges(a))
    priority[b] = Inf
  }
  if (verbose) cat('\n')
  
  if (! missing(save.path))
    save(cg, file = paste0(save.path,
                           'contractible_graph_',
                           num.components, '_',
                           Sys.Date(),
                           '.RData'))

  if (return.graph) return(cg)
  cg$get_vertex_components()
}

"
edges = list()
edges$edge.mat = t(matrix(c(1, 2,
                            1, 3,
                            2, 1,
                            2, 4,
                            3, 1,
                            3, 4,
                            4, 2,
                            4, 3),
                          nrow = 2))
edges$energy.vec = c(0.7, 0.8, 0.7, 0.5, 0.8, 0.9, 0.5, 0.9)
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
