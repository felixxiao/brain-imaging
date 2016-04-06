# compute criteria
"
source('criteria/header_criteria.R')
is_connected(edges_nz$edge.mat, part$bf.ec)
criterion.adjacent_pairwise_ecor(edges_nz, part$bf.ec)$total.mean
criterion.cut_weight(edges_nz, part$bf, 'ratio')
criterion.balance(part$bf.ec)
criterion.jaggedness(edges_nz, part$bf.ec)
"
