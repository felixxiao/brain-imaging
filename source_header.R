#load all the functions in all the source files in this entire git repo

#note: some functions are hidden. to 
# show them, use "ls(all.names = TRUE)"

script.dir = dirname(parent.frame(2)$ofile)
script.dir = normalizePath(script.dir)

#add folders to this vector
folder.vec = c(
 "common_func",
 "energy_parcellation",
 "parcel_criteria",
 "preprocessing_script", 
 "plotter",
 "union_find")
setwd(script.dir)

for(i in 1:length(folder.vec)){
  setwd(folder.vec[i])

  #find the one file with the keyword "header"
  source(grep("header", dir(), value = T)[1])

  setwd(script.dir)
}

#remove the extra variables
rm(list = c("script.dir", "folder.vec", "i"))

