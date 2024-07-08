## Function to pull out the gene target of a drug that has the max correlation.
PredTarget <- function(Sim.GES.DRS,D.M){
    ## Make sure the user use the same meta data for the list of Sim.GES.DRS
    Drug.Id.i <- match(names(Sim.GES.DRS), D.M[,1])
    D.M.f <- D.M[Drug.Id.i,]
    ## based on the name of the drug.
    Splt.targets <- str_split(as.character(D.M.f[,3]), ", ")
    F.V <- nrow(Sim.GES.DRS[[1]])
    corrMatPval <- vapply(Sim.GES.DRS, function(x) x[,1],numeric(F.V))
    corrMat <- vapply(Sim.GES.DRS, function(x) x[,2],numeric(F.V))
    corrMatFDR <- vapply(Sim.GES.DRS, function(x) x[,3],numeric(F.V))
    Stat.Drug <- bplapply(seq_len(length(Splt.targets)),function(x){
        ret_cor <- corrMat[
        match(Splt.targets[[x]],
        rownames(corrMat)), x]
        names(ret_cor) <- tryCatch(Splt.targets[[x]],error=function(e){NA})
        ret_cor
        })
    ## if there is one drug.
    if ( nrow(D.M.f)==1){
        ## Catch potential error of returning -Inf
        Target.Max.cor <- tryCatch(max(Stat.Drug,na.rm=TRUE ),
            error=function(e){NA})
        Target.Max.Name <- names(Stat.Drug)[which.max(Stat.Drug)]
        ## drug Id will be the same as the metadata.
        Drug.id <- names(Sim.GES.DRS)
        names(Target.Max.cor) <- as.character(D.M.f[,2])
    }else{
        names(Stat.Drug) <- as.character(D.M.f[,2])
        Target.Max.cor <- vapply(Stat.Drug,
            function(x) max(x, na.rm=TRUE),numeric(1))
        Target.Max.Name <- vapply(Stat.Drug,
            function(x) names(which.max(x)),character(1))
        Drug.id <- D.M.f[match(names(Target.Max.cor), D.M.f$name),1]
    }
    Target.Max.Name[vapply(Target.Max.Name, length,FUN.VALUE =numeric(1))==0] <- NA
    ## Obtain drug ID and make sure that we map the correct one
    Target.Pred <- data.frame(
        DrugID = Drug.id,
        drugName=names(Target.Max.cor),
        MaxTargetName=unlist(Target.Max.Name),
        Maxcorr=Target.Max.cor)
    ## Pull out p val and FDR.
    Target.Pred$KnownTargetCorrP <- vapply(seq_len(nrow(Target.Pred)), function(x)
        tryCatch(corrMatPval[Target.Pred[x,3], Target.Pred[x,1]],
        error=function(e){NA}),FUN.VALUE = numeric(1))
    # Known Target Significance - FDR corrected
    Target.Pred$KnownTargetCorrFDR <- vapply(seq_len(nrow(Target.Pred)), function(x)
        tryCatch(corrMatFDR[Target.Pred[x,3], Target.Pred[x,1]],
        error=function(e){NA}),FUN.VALUE = numeric(1))
    Target.Pred
}


