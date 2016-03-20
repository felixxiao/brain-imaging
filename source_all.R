script.dir = dirname(sys.frame(1)$ofile)
script.dir = normalizePath(script.dir)

#add folders to this vector
folder.vec = c(
 "energy_parcellation",
 "criteria",
 "partition_edgecontract", 
 "partition_opt",
 "plotter",
 "preprocess"
)

for(i in folder.vec){
  #find the one file with the keyword "header"
  name = paste0(i, '/', grep("^header_", dir(i), value = T))
  if(length(name) != 1) stop(paste0("Too many header files in", paste0))

  source(name)
}

#remove the extra variables
rm(list = c("script.dir", "folder.vec", "i", "name"))
