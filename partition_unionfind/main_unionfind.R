setwd('~/GitHub/union_find')
source('header_unionfind.R')

load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

part1 = partition.addedge.uf(edges, max.size = 15000, min.size = 1000, PATH_SAVE)
#load(paste0(PATH_DATA, 'edge_connect_1000_15000_2016-01-19.RData'))
table(part1)

part2 = partition.addedge.uf(edges, max.size = 10000, min.size = 1000, PATH_SAVE)
table(part2)

part3 = partition.addedge.uf(edges, max.size = 7500, min.size = 1000, PATH_SAVE)
table(part3)


setwd('~/GitHub/parcel_criteria')
source('header_criteria.R')
load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
#load(paste0(PATH_DATA, 'template_2015-12-07.RData'))

crit1 = list()
crit2 = list()
crit3 = list()

crit1$within = criterion.within_pairwise_ecor(dat$dat, part1)
crit1$between = criterion.between_pairwise_ecor(dat$dat, part1)

crit2$within = criterion.within_pairwise_ecor(dat$dat, part2)
crit2$between = criterion.between_pairwise_ecor(dat$dat, part2)

crit3$within = criterion.within_pairwise_ecor(dat$dat, part3)

# compute confidence interval
object_size(crit1)
