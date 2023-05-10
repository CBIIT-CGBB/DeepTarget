## This is used to pull out the gene that has the max correlation with the drug of interest.

PredMaxSim <- function( Sim.GES.DRS=Sim.GES.DRS, D.M = Drug.Metadata){
    Drug.Id.i <- match(names(Sim.GES.DRS), D.M[,1])
    D.M.f <- D.M[Drug.Id.i,]
    ## Based on the name of the drug:
    Splt.targets = str_split(as.character(D.M.f[,3]), ", ")
    ## This is extracts P values and correlation values, and turns them into a matrix where the genes are rows and the drugs are columns.
    corrMatPval=sapply(Sim.GES.DRS, function(x) x[,1])
    corrMat=sapply(Sim.GES.DRS, function(x) x[,2])
    ## Calculate false discovery rate (FDR).
    corrMatFDR=sapply(Sim.GES.DRS, function(x) x[,3])
##############
    ## If there is one drug:
    if ( nrow(D.M.f)==1){
        BestTargetGene=rownames(corrMat)[which.max(corrMat[,1])]
        BestTargetCorr = max(corrMat[,1],na.rm = T)
        Drug.id <- names(Sim.GES.DRS)
    }else{
        ## 2 is the column.
        BestTargetGene=apply(corrMat,
                             2, function(x) rownames(corrMat)[which.max(x)] )
        BestTargetCorr=apply(corrMat,2, function(x) max(x, na.rm = T) )
        Drug.id = names(BestTargetGene)
    }


Pred.sim=data.frame(
    DrugID = Drug.id,
    BestTargetGene=unlist(BestTargetGene),
    BestTargetCorr=BestTargetCorr)

# Map to get the P value from gene that is the best target.
Pred.sim$BestTargetCorrP = sapply(1:nrow(Pred.sim), function(x)
    errHandle(corrMatPval[Pred.sim[x,'BestTargetGene'], Pred.sim[x,1]]) )
# Best Hit Significance - FDR corrected
Pred.sim$BestTargetCorrFDR = sapply(1:nrow(Pred.sim), function(x)
    errHandle(corrMatFDR[Pred.sim[x,'BestTargetGene'], Pred.sim[x,1]]) )

Pred.sim

}

