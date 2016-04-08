files = grep('-partitions.RData', dir('ABIDE'), value = T)
files = paste0('ABIDE/', files)

brains = grep('-matrix.RData', dir('ABIDE'), value = T)
brains = sapply(strsplit(brains, '[[:punct:]]'), function(x) x[1])
brains = paste0('ABIDE/', brains)

#### Generate random parcellations of each brain
source('partition_edgecontract/header_contraction.R')

part.rand = lapply(brains, function(brain) {
  load(paste0(brain, '-matrix.RData'))
  load(paste0(brain, '-edges.RData'))
  edges_nz = preprocess.all(dat$dat, edges)$edges_nz
  edges_nz$energy.vec = sample(edges_nz$energy.vec)
  partition.contractedge.general(116, 6, 4, edges_nz)$partitions[[1]]
})
names(part.rand) = brains

#### Compare different parcellations of the same brain

library(mclust)
library(gplots)

ari = sapply(brains, function(brain) {
  load(paste0(brain, '-partitions.RData'))
  part[['Shuffled GenEC (6,4)']] = part.rand[[brain]]
  sapply(part, function(p) sapply(part, adjustedRandIndex, x = p))
}, simplify = 'array')
dimnames(ari)[[1]] = dimnames(ari)[[2]] =
  c('GenEC (3,1)', 'GenEC (6,1)', 'GenEC (6,4)', 'Spectral',
    'Spectral GenEC (6,4)', 'SymBMF', 'SymBMF GenEC (6,4)', 'Shuffled GenEC (6,4)')
ari.mean = apply(ari, c(1,2), mean)

colors = RColorBrewer::brewer.pal(9, name = 'Blues')
png('writeup/figs/8_ari_same.png', width = 1000, height = 800, res = 150)
par(oma = c(5, 4, 1, 9))
heatmap.2(ari.mean, dendrogram = 'none', Rowv = F, Colv = F, scale = 'none',
          cellnote = round(ari.mean, 3), trace = 'none', col = colors,
          key = F, keysize = 0, notecol = 'red', srtCol = 45)
dev.off()

#### Compare across different brains
source('energy_parcellation/header_energy.R')

maps = lapply(brains, function(brain) {
  load(paste0(brain, '-matrix.RData'))
  load(paste0(brain, '-edges.RData'))
  preprocess.all(dat$dat, edges)$map
})

part = lapply(brains, function(brain) {
  load(paste0(brain, '-partitions.RData'))
  part$ec.6.4
})

m = length(part)

ari = matrix(NA, nrow = m, ncol = m)
for (i in 2:m)
  for (j in 1:(i-1))
    ari[i,j] = ari[j,i] = adjustedRandIndex(part[[i]][maps[[i]] %in% maps[[j]]],
                                            part[[j]][maps[[j]] %in% maps[[i]]])
rownames(ari) = colnames(ari) = sapply(strsplit(brains, '/'), function(x) x[2])
colors = RColorBrewer::brewer.pal(9, name = 'Blues')
png('writeup/figs/8_ari_diff.png', width = 1500, height = 1200, res = 200)
par(oma = c(3, 2, 1, 4))
heatmap.2(ari, dendrogram = 'none', Rowv = F, Colv = F, scale = 'none',
          cellnote = round(ari, 3), trace = 'none', col = colors,
          key = F, keysize = 0, notecol = 'red', notecex = 0.8, srtCol = 45)
dev.off()

mean(ari[1:6,7:12])
(mean(ari[1:6,1:6], na.rm = T) + mean(ari[7:12,7:12], na.rm = T)) / 2
