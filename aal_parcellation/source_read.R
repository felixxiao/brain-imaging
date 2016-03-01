#reads in a filename (the AAL file) and then output
# the corresponding partitions based on the AAL
#the inputs of adj.list and mask are from the 
# template
read.aal <- function(filename, adj.list, mask){
  edge.mat = convert.adjList2edgeMat(adj.list)
  
  #read in the data
  aal = readNIfTI(filename)
  aal = aal@.Data

  #find out the "zero.voxels" which are ones in
  # mni but not in aal
  idx.aal = which(aal != 0)
  idx.aal = idx.aal[idx.aal %in% mask] 
  
  #form the base.line graph
  zero.voxel = which(!mask %in% idx.aal)
  n = max(edge.mat)
  g = .construct.graphBase(zero.voxel, edge.mat, n)
  
  #add in all the parcellations in aal
  aal.val = as.numeric(aal)[idx.aal]
  uniq.aal = unique(aal.val)

  for(i in 1:length(uniq.aal)){
    idx.parcel = which(aal.val == uniq.aal[i])
    voxels.inParcel = idx.aal[idx.parcel]

    #find the relevant edges in edge.mat
    leftRight.bool = apply(edge.mat, 2, function(x){
      which(x %in% voxels.inParcel)
    })
    total.bool = apply(leftRight.bool, 1, all)
    edges.inParcel = which(total.bool == TRUE)
   
    #add the selected edges into the graph g
    g = add.edges(g, t(edge.mat[edges.inParcel,]))
  }

  components(g)$no
}
