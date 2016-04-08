source('energy_parcellation/header_energy.R')
#brain = commandArgs()
load(paste0(brain, '-partitions.RData'))

# get dat and edges for brain
load(paste0(brain, '-matrix.RData'))
load(paste0(brain, '-edges.RData'))
dat = dat$dat
N = ncol(dat)

# form dat_nz and edges_nz
out = preprocess.all(dat, edges)
dat_nz = out$dat_nz
edges_nz = out$edges_nz
map_nz = out$map_nz
rm(out)

source('criteria/header_criteria.R')

connect = sapply(part, is_connected, edge.mat = edges_nz$edge.mat)
adjcent = sapply(part, function(p) criterion.adjacent_pairwise_ecor(edges_nz, p)$total.mean)
ratioct = sapply(part, criterion.cut_weight, edges = edges_nz, type = 'ratio')
balance = sapply(part, criterion.balance)
jaggedn = sapply(part, criterion.jaggedness, edges = edges_nz)
boundry = sapply(part, criterion.multi_boundary_ecor, dat = dat_nz)

rm(dat_nz, edges_nz, map_nz, N, dat, edges)
save(list = ls(), file = paste0(brain, '-criteria.R'))

"
crit = cbind(connect, adjcent, boundry, ratioct, balance, jaggedn)
colnames(crit) = c('Connected', 'Adjacent', 'Boundary', 'Ratio-Cut', 'Balance', 'Jaggedness')
rownames(crit) = c('GenEC (3,1)', 'GenEC (6,1)', 'GenEC (6,4)', 'Spectral',
                   'Spectral GenEC (6,4)', 'SymBMF', 'SymBMF GenEC (6,4)')

library(gplots)
crit.text = apply(crit, c(1,2), signif, digits = 4)
crit.text[,1] = plyr::mapvalues(crit.text[,1], c(0, 1), c('No', 'Yes'))
colors = RColorBrewer::brewer.pal(9, name = 'Blues')
png('crit.png', width = 1200, height = 800, res = 200)
par(oma = c(3, 1, 1, 1))
heatmap.2(crit, dendrogram = 'none', Rowv = F, scale = 'column',
          cellnote = crit.text, trace = 'none', col = colors,
          key = F, keysize = 0, notecol = 'red')
dev.off()
"