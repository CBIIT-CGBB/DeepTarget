## This is the pathway analysis based on the information from the drug such as MOA pathways and genes targeted.
DoPWY <- function(Sim.GES.DRS=Sim.GES.DRS, D.M = Drug.Metadata){
  PwyTargeted=split(D.M[,3], D.M[,5])
  PwyTargeted.unl=sapply(PwyTargeted, unlist)
  PwyTargeted.c=lapply(PwyTargeted.unl, function(y){
    unique(unlist(sapply(as.character(y), function(x)
      unlist(strsplit(x, ', ')) )))
  } )

  PwyTargeted.c=sapply( PwyTargeted.c, na.omit)
  corrMat=sapply(Sim.GES.DRS, function(x) x[,2])
  #drugVScrispr_corr_Strength=sapply(infunc_drugVScrispr_corr_features_list, function(x) x[,2])
  DoPWYEnr=mclapply(1:ncol(corrMat), function(x) {
    corrMat.u=unlist(corrMat[,x])
    names(corrMat.u)  = rownames(corrMat)
    corrMat.u.O=sort(corrMat.u, decreasing = T)
    gseaEnr=fgsea(pathways = PwyTargeted.c,
                          stats = corrMat.u.O,
                          minSize=1,
                          maxSize=100)
    gseaEnr
  }, mc.cores = 2)
  names(DoPWYEnr)=colnames( corrMat)
  DoPWYEnr
}
