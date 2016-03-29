setwd("~/felix_senior_thesis_2015-16/code.git")

source("source_all.R")

setwd("~/Felix_Thesis_2016/ABIDE")
load("../template_2016-03-26.RData")

files = dir()
files = files[grep("ABIDE.*-matrix.*.RData", files)]

edge.mat = convert.adjList2edgeMat(template$neighbor.list)

for(i in 1:length(files)){
  load(files[i])

  assert_that(length(template$neighbor.list) == ncol(dat$dat))

  #get the subject ID
  subjID = strsplit(files[i], "-")[[1]][2]

  edges = compute.edgeWeights(dat$dat, template$neighbor.list, edge.mat = edge.mat, save = F,
    parallel = T)
  save(edges, file = paste0("ABIDE-", subjID, "-edgeWeights_", DATE, ".RData"))
}
