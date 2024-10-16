## Function to execute the pathway analysis for a set of targets.
DoPWY <- function(Sim.GES.DRS,D.M){
    PwyTargeted <- split(D.M[,3], D.M[,5])
    PwyTargeted.unl <- lapply(PwyTargeted, unlist)
    PwyTargeted.c <- lapply(PwyTargeted.unl, function(y){
        unique(unlist(lapply(as.character(y), function(x)
            unlist(strsplit(x, ', ')))))
    })
    PwyTargeted.c <- bplapply( PwyTargeted.c, na.omit)
    F.V <- nrow(Sim.GES.DRS[[1]])
    corrMat <- vapply(Sim.GES.DRS, function(x) x[,2],numeric(F.V))
    DoPWYEnr <- bplapply(seq_len(ncol(corrMat)), function(x) {
        corrMat.ul <- unlist(corrMat[,x])
        names(corrMat.ul) <- rownames(corrMat)
        corrMat.ul.s <- sort(corrMat.ul, decreasing = TRUE)
        corrMat.ul.s <- as.matrix(corrMat.ul.s)
        ## remove the duplicated gene if having
        corrMat.ul.s.u <- corrMat.ul.s[!duplicated(row.names(corrMat.ul.s)),]
        gseaEnr <- fgsea(
            pathways = PwyTargeted.c,stats = corrMat.ul.s.u,minSize=1,maxSize=100)
        gseaEnr})
    names(DoPWYEnr) <- colnames(corrMat)
    DoPWYEnr
}

