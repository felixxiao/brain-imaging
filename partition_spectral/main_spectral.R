source('energy_parcellation/header_energy.R')

load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
load(paste0(PATH_DATA, 'template_2015-12-07.RData'))

N = ncol(dat$dat)

rm_zero = preprocess.remove_zero_vertices(dat$dat)
dat = rm_zero$data
map_nonzero = rm_zero$map
rm(rm_zero)

edge.mat = convert.adjList2edgeMat(template$neighbor.list, duplicates = F)
edge.mat = preprocess.map_vertices(edge.mat, map_nonzero, inverse = TRUE)
edge.mat = edge.mat[!(is.na(edge.mat[,1]) | is.na(edge.mat[,2])),]
#edges = compute.edgeWeights(dat, NULL, edge.mat = edge.mat, verbose = T)
text = 'ABIDE 50002 with zero columns removed'
save(edges = edges, details = text, N = N, map_nonzero = map_nonzero,
     file = paste0(PATH_DATA, 'edges_2016-01-29.RData'))


load(paste0(PATH_DATA, 'edges_2016-01-29.RData'))

source('partition_opt/header_opt.R')
part.rc = partition.spectral.recursive(edges, max.depth = 7)
part.mt = partition.spectral.multiway(edges, 128)

map.back = function(part, map, N)
{
  part.map = factor(rep(NA, times = N), levels = 1:(nlevels(part) + 1))
  part.map[map] = part
  part[is.na(part)] = nlevels(part) + 1
  part.map
}

part.rc_map = map.back(part.rc, map_nonzero, N)
part.mt_map = map.back(part.mt, map_nonzero, N)

source('criteria/header_criteria.R')
load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))
adj = criterion.adjacent_pairwise_ecor(edges, part.rc_map)
bou = criterion.boundary_pairwise_ecor(edges, part.rc_map)

source('plotter/header_plotter.R')
MNI = load.nifti(paste0(PATH_DATA, 'MNI152_T1_2mm_brain_symmetric.nii.gz'))
plot.partition2D(part.mt_map, MNI, template$mask,
                 filename = paste0('spectral_multi_', nlevels(part.mt_map)),
                 view = VIEWS[3])
