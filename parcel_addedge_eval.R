source('union_find/header_unionfind.R')
source('parcel_criteria/header_criteria.R')

load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

# generate edge list with shuffled energy cov values
edges.rand = edges
edges.rand$energy.vec = sample(edges$energy.vec, length(edges$energy.vec))
# parcellate random edges
parcel.rand = partition.addedge.uf(edges.rand, max.size = 10000, min.size = 1000)

save(list = ls(), paste0(PATH_SAVE, 'parcel_addedge_eval_', DATE, '.RData'))

# max and min component size parameters for add edge algorithm
size = t(matrix(c(15000, 1000,
                  10000, 1000,
                   7500, 1000,
                   5000, 1000,
                  10000, 1500,
                  10000,  750,
                  10000,  500),
                ncol = 2))
colnames(size) = c('max', 'min')

# parcellate with add edge algorithm with above size parameters
parcel = list()
for (i in 1:nrow(size))
{
  parcel[[i]] = partition.addedge.uf(edges, size[i,'max'], size[i,'min'])
  save(list = ls(), paste0(PATH_SAVE, 'parcel_addedge_eval_', DATE, '.RData'))
}

# now evaluate the parcellations
criteria = list()
for (i in 1:nrow(size))
{
  criteria[[i]] = list()
  criteria[[i]]$within   = criterion.within_pairwise_ecor(dat$dat, parcel[[i]], subsample = 100)
  criteria[[i]]$between  = criterion.between_pairwise_ecor(dat$dat, parcel[[i]], subsample = 100)
  criteria[[i]]$adjacent = criterion.adjacent_pairwise_ecor(edges, parcel[[i]])
  criteria[[i]]$boundary = criterion.boundary_pairwise_ecor(edges, parcel[[i]])
  save(list = ls(), paste0(PATH_SAVE, 'parcel_addedge_eval_', DATE, '.RData'))
}

criteria.rand$within   = criterion.within_pairwise_ecor(dat$dat, parcel.rand, subsample = 4000)
criteria.rand$between  = criterion.between_pairwise_ecor(dat$dat, parcel.rand, subsample = 4000)
criteria.rand$adjacent = criterion.adjacent_pairwise_ecor(edges, parcel.rand)
criteria.rand$boundary = criterion.boundary_pairwise_ecor(edges, parcel.rand)
save(list = ls(), paste0(PATH_SAVE, 'parcel_addedge_eval_', DATE, '.RData'))

#load(paste0(PATH_DATA, 'parcel_addedge_eval_2016-01-23.RData'))

sink('writeup/3_constrained_addedge.txt')
for (i in 1:length(criteria))
{
  cat(sprintf('%5d', size[i, c('min', 'max')]), sep = ' & ')
  cat(' & ')
  crit = sapply(c('within', 'adjacent', 'between', 'boundary'),
                function(x) criteria[[i]][[x]]$total.mean)
  crit = sprintf('%.3f', crit)
  cat(crit, sep = ' & ')
  cat('\\\\ \n')
}
cat('\\hline\n')
cat('& Random &')
cat(sprintf('%.3f', sapply(c('within', 'adjacent', 'between', 'boundary'),
                           function(x) criteria.rand[[x]]$total.mean)),
    sep = ' & ')
cat('\n')
sink()

load(paste0(PATH_DATA, 'template_2015-12-07.RData'))
source('plotter/header_plotter.R')
MNI = load.nifti(paste0(PATH_DATA, 'MNI152_T1_2mm_brain_symmetric.nii.gz'))

plot.partition2D(parcel[[3]], MNI, template$mask,
                 'sc_addedge_1000_7500', view = VIEWS[3])
