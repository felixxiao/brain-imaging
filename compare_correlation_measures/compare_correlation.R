source('energy_parcellation/header_energy.R')

load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))
load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))

edges$pearson = compute.edgeWeights(dat$dat, NULL, cor, edge.mat = edges$edge.mat,
                                    verbose = T, save = F)$energy.vec

all((edges$energy.vec == 0) == is.na(edges$pearson))

idx.rm = (edges$energy.vec == 0)

edges$energy.vec[idx.rm] = Inf
edges$pearson[idx.rm]    = Inf

energy.rank  = rank(edges$energy.vec)
pearson.rank = rank(edges$pearson)
r2.rank = rank(edges$pearson ^ 2)

energy.rank[idx.rm]  = 0
pearson.rank[idx.rm] = 0
r2.rank[idx.rm ] = 0

hist(energy.rank - pearson.rank)
order(energy.rank - pearson.rank, decreasing = F)[1:10]
order(energy.rank - r2.rank, decreasing = F)[1:10]

E = c(363027, 873267, 244037, 507420, 265147) # ecor > r2 [rankwise] b/c non-linear
E = c(E, 382613, 696479) # ecor < r2 [rankwise] b/c outliers & linearity assumption

png('writeup/1_nonlinear_ABIDE_50002.png', width = 450, height = 500)
par(mfrow = c(2,2))
for (e in E[1:4])
{
  i = edges$edge.mat[e,1]; j = edges$edge.mat[e,2]
  plot(dat$dat[,i], dat$dat[,j],
       xlab = paste('Voxel', i),
       ylab = paste('Voxel', j))
}
dev.off()