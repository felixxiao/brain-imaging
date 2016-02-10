# TODO: spectral partition into k-partitions

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
load(paste0(PATH_DATA, 'edges_2016-01-29.RData'))

source('partition_spectral/source_spectral.R')
part.sp = partition.spectral(edges)
assert_that(length(part.sp) == length(map_nonzero))

part.sp_map = factor(rep(NA, times = N), levels = c(1,2,3))
part.sp_map[map_nonzero] = part.sp
part.sp_map[is.na(part.sp_map)] = 3

source('parcel_criteria/header_criteria.R')
load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))
adj = criterion.adjacent_pairwise_ecor(edges, part.sp_map)
bou = criterion.boundary_pairwise_ecor(edges, part.sp_map)

source('plotter/header_plotter.R')
MNI = load.nifti(paste0(PATH_DATA, 'MNI152_T1_2mm_brain_symmetric.nii.gz'))
plot.partition2D(part.sp_map, MNI, template$mask,
                 filename = 'spectral', view = VIEWS[1])

