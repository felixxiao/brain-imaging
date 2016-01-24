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
    # - edges has duplicates
    # - in edges$edge.mat, vertices are numbered 1:n
    initialize = function(edges)
    {
      edge.mat = edges$edge.mat
      self$n = max(edge.mat)
      energy.vec = edges$energy.vec
      assert_that(length(energy.vec) == nrow(edge.mat))
      
      n = max(edge.mat)
      self$vertices = as.list(1:n)
      names(self$vertices) = as.character(1:n)
      self$edges = lapply(1:n, function(x) list())
      names(self$edges) = as.character(1:n)
      
      for (e in 1:nrow(edge.mat))
      {
        a = as.character(edge.mat[e,1])
        b = as.character(edge.mat[e,2])
        self$edges[[a]][[b]] = energy.vec[e]
        if (! is.null(self$edges[[b]][[a]]))
          assert_that(self$edges[[b]][[a]] == self$edges[[a]][[b]])
      }
    },
    
    get_edges = function(comp)
    {
      c = as.character(comp)
      weights = sapply(self$edges[[c]], mean)
      names(weights) = names(self$edges[[c]])
      weights
    },
    
    get_size = function(comp)
    {
      c = as.character(comp)
      length(self$vertices[[c]])
    },
    
    contract_components = function(comp1, comp2)
    {
      a = as.character(comp1)
      b = as.character(comp2)
      
      if (is.null(self$edges[[a]][[b]]) | is.null(self$edges[[b]][[a]]))
      {
        assert_that(is.null(self$edges[[b]][[a]]) & is.null(self$edges[[a]][[b]]))
        warning(comp1, ' and ', comp2, ' are not connected')
        return()
      }
      
      assert_that(length(intersect(self$vertices[[a]], self$vertices[[b]])) == 0)
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
    
    get_components = function(as.factor = T)
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

partition.contractedge = function(edges, num.components)
{
  cg = ContractibleGraph$new(edges)
  n = max(edges$edge.mat)
  assert_that(num.components < n)
  
  components = rep(NA, times = n)
  for (c in 1:n)
    components[c] = 1 - max(cg$get_edges(c))
  names(components) = as.character(1:n)
  
  for (i in 1:(n - num.components))
  {
    a = names(which.min(components))
    b = names(which.max(cg$get_edges(a)))
    c = cg$contract_components(a, b)
    assert_that(c == a)
    components[a] = cg$get_size(a) - max(cg$get_edges(a))
    components[b] = Inf
  }
  
  cg$get_components()
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
