#note: some functions are hidden. to 
# show them, use "ls(all.names = TRUE)"

script.dir = dirname(parent.frame(2)$ofile)
script.dir = normalizePath(script.dir)

#add folders to this vector
folder.vec = c(
 "preprocessing_script", 
 "energy_parcellation")
setwd(script.dir)

for(i in 1:length(folder.vec)){
  setwd(folder.vec[i])
  source("source_header.R")
  setwd(script.dir)
}

