source('partition_unionfind/header_unionfind.R')

load(paste0(PATH_DATA, 'results/edges_ABIDE50002_nz.RData'))

uf1 = partition.addedge.uf(edges_nz, max.size = 5000, min.size = 100,
                           min.components = 116, return.uf = T)
part.uf = uf1$get_parcels()
plot(as.vector(table(part.uf)),
     criterion.adjacent_pairwise_ecor(edges_nz, part.uf)$mean)
as.vector(table(part.uf))
i = which(part.uf == 1)

edges_nz = preprocess.validate.edges(edges_nz)
get.edges(edges_nz, i[2], NULL)

get.edges = function(edges, v, w, summary.func = NULL)
{
  .validate.edges(edges)
  if (! is.null(w))
  {
    assert_that(v < w)
    e = (edges$edge.mat[,1] == v & edges$edge.mat[,2] == w)
  }
  else
  {
    e = (edges$edge.mat[,1] == v | edges$edge.mat[,2] == v)
  }
  if (is.null(summary.func))
    return(edges$energy.vec[e])
  summary.func(edges$energy.vec[e])
}

