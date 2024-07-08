## Function to predict the maximum similarity for gene targets.
PredMaxSim <- function(Sim.GES.DRS,D.M){
    Drug.Id.i <- match(names(Sim.GES.DRS), D.M[,1])
    D.M.f <- D.M[Drug.Id.i,]
    ## based on the name of the drug.
    Splt.targets <- str_split(as.character(D.M.f[,3]), ", ")
    F.V <- nrow(Sim.GES.DRS[[1]])
    corrMatPval <- vapply(Sim.GES.DRS, function(x) x[,1],numeric(F.V))
    corrMat <- vapply(Sim.GES.DRS, function(x) x[,2],numeric(F.V))
    corrMatFDR <- vapply(Sim.GES.DRS, function(x) x[,3],numeric(F.V))
    ## if there is one drug.
    if ( nrow(D.M.f)==1){
        BestTargetGene <- rownames(corrMat)[which.max(corrMat[,1])]
        BestTargetCorr <-  max(corrMat[,1],na.rm = TRUE)
        Drug.id <- names(Sim.GES.DRS)
    }else{
        ## 2 is the column.
        BestTargetGene <- apply(
            corrMat,2, function(x) rownames(corrMat)[which.max(x)] )
        BestTargetCorr <- apply(
            corrMat,2, function(x) max(x, na.rm = TRUE) )
        Drug.id <- names(BestTargetGene)
    }
    Pred.sim <- data.frame(
        DrugID = Drug.id, BestTargetGene=unlist(BestTargetGene),
        BestTargetCorr=BestTargetCorr)
    ## map to get the P val from best target gene
    Pred.sim$BestTargetCorrP <- vapply(seq_len(nrow(Pred.sim)), function(x)
    tryCatch(corrMatPval[Pred.sim[x,'BestTargetGene'], Pred.sim[x,1]], 
        error=function(e){NA}),numeric(1))
    ## Best Hit Significance - FDR corrected
    Pred.sim$BestTargetCorrFDR <- vapply(seq_len(nrow(Pred.sim)),function(x) 
        tryCatch(corrMatFDR[Pred.sim[x,'BestTargetGene'], Pred.sim[x,1]],
            error=function(e){NA}), numeric(1))
    Pred.sim
}

