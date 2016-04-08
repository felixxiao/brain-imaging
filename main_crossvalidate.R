source('energy_parcellation/header_energy.R')
source('criteria/header_criteria.R')

cross_validate = function(part_nz, dat, edges, edges_test)
{
  map_nz = preprocess.all(dat, edges)$map_nz
  N = max(edges$edge.mat)
  part = preprocess.assign_unmapped_vertices(map_nz, edges_test, N, part_nz)
  criterion.adjacent_pairwise_ecor(edges_test, part)$total.mean
}

brains = sapply(strsplit(grep('-matrix.RData', dir('ABIDE'), value = T), '-'),
                function(x) x[1])

cv = matrix(NA, nrow = length(brains), ncol = length(brains), 
            dimnames = list(brains, brains))

for (brain.part in brains)
  for (brain.test in brains)
  {
    load(paste0('ABIDE/', brain.part, '-partitions.RData'))
    load(paste0('ABIDE/', brain.test, '-edges.RData'))
    edges_test = edges
    load(paste0('ABIDE/', brain.part, '-edges.RData'))
    load(paste0('ABIDE/', brain.part, '-matrix.RData'))
    
    cv[brain.part, brain.test] = cross_validate(part$ec.3.1, dat$dat,
                                                edges, edges_test) -
                                 mean(edges_test$energy.vec)
  }

library(gplots)
colors = RColorBrewer::brewer.pal(11, name = 'RdBu')
png('writeup/figs/8_cv_ec.png', width = 1500, height = 1200, res = 200)
par(oma = c(2, 1.5, 1, 2))
heatmap.2(cv, dendrogram = 'none', Rowv = F, Colv = F, scale = 'none',
          cellnote = round(cv, 3), trace = 'none', col = colors,
          key = F, keysize = 0, notecol = 'red', notecex = 0.8, srtCol = 45,
          xlab = 'Parcellation Brain', ylab = 'Validation Brain',
          offsetRow = 0, offsetCol = 0)
dev.off()

cv.no_diag = cv
diag(cv.no_diag) = NA
hist(cv.no_diag)

