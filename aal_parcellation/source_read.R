#reads in a filename (the AAL file) and then output
# the corresponding partitions based on the AAL
#the inputs of adj.list and mask are from the 
# template
#enforce.connectivity is true if you want to "lazily" enforce
# connectivity within each parcel. note that if false, you'll
# get more than 116 parcels since some parcels aren't "technically"
# connected. 
read.aal <- function(filename, adj.list, mask, enforce.connectivity = TRUE,
 verbose = FALSE){
  edge.mat = convert.adjList2edgeMat(adj.list)
  
  ##REPLACE WITH LOAD.NIFTI
  #read in the data
  aal = readNIfTI(filename)
  aal = aal@.Data

  #find out the "zero.voxels" which are ones in
  # mni but not in aal
  idx.aal = which(aal != 0)
  idx.aal = idx.aal[idx.aal %in% mask] 
  idx.aalconvert = mapvalues(idx.aal, mask, 
   1:max(edge.mat), warn_missing = F)  
  #note: idx.aal are voxel locations in 3D, but idx.aalconvert
  # are values that go from 1:max(edge.mat)

 #form the base.line graph
  zero.voxels = which(!mask %in% idx.aal)
  n = max(edge.mat)
  g = .construct.graphBase(zero.voxels, edge.mat, n)
  
  #add in all the parcellations in aal
  aal.val = as.numeric(aal)[idx.aal]
  uniq.aal = unique(aal.val)

  old.comp = components(g)$no

  for(i in 1:length(uniq.aal)){
    idx.parcel = which(aal.val == uniq.aal[i])
    voxels.inParcel = idx.aalconvert[idx.parcel]


    if(enforce.connectivity){
      #simply add each edge to each other in a chain graph
      edges.toadd = matrix(0, nrow = length(voxels.inParcel) - 1, ncol = 2)
      edges.toadd[,1] = voxels.inParcel[1:(length(voxels.inParcel)-1)]
      edges.toadd[,2] = voxels.inParcel[2:length(voxels.inParcel)]

    } else {
      #find the relevant edges in edge.mat
      leftRight.bool = alply(edge.mat, 2, function(x){
        which(x %in% voxels.inParcel)
      })
    
      edges.inParcel = intersect(leftRight.bool[[1]], leftRight.bool[[2]])
      edges.toadd = edge.mat[edges.inParcel,]
    }  

    #add the selected edges into the graph g
    g = add.edges(g, t(edges.toadd))
  
    if(verbose){
      #check if the components went down the right amount
      new.comp = components(g)$no
      if(old.comp - new.comp != length(idx.parcel)-1){
        print(paste0("PARCEL: ", uniq.aal[i], " for i = ", i, " with supposed decrease in ",
         length(idx.parcel)-1, " but actual drop in ", old.comp - new.comp))
      }
      old.comp = new.comp

      if(i %% floor(length(uniq.aal)/10) == 0) cat('*')
    }
  }

  as.factor(components(g)$membership)
}
