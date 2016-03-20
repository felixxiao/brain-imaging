setwd('~/../OneDrive/data-night1')
source('lib.R')
setwd('~/GitHub')
source('partition_nmf/source_nmf.R')
source('partition_opt/header_opt.R')

k = 15
n = sample(5:65, size = k, replace = T)
A = generate.groups(n, 0.6, 0.1, 20, 5)
S = A$S
A = A$A
edges = compute.edges(A)

# Partitioning by SymNMF
symnmf = partition.sym_nmf(A, k, 1, .04, iter = 200)
plot(symnmf$obj, type = 'l')
sum(apply(symnmf$H, 1, which.max) == apply(symnmf$W, 1, which.max)) / nrow(symnmf$H)
S1 = split(1:nrow(symnmf$H), apply(symnmf$H + symnmf$W, 1, which.max))

# Partitioning by Spectral
S2 = partition.spectral.multiway(edges, k)
S2 = split(1:length(S2), S2)

#sapply(S2, length)
#plot.network(A, S2)
#empirical.p(A, S1)

sort(n)
sort(unname(sapply(S2, length)))

cut.weight(A, S)
cut.weight(A, S1)
cut.weight(A, S2)

rand.index(S, S1, sum(n))
rand.index(S, S2, sum(n))
