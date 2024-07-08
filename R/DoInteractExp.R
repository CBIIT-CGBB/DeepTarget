## Function to calculate the WT target interaction based on a cut-off value.
DoInteractExp <- function(Predtargets,Exp,DRS,GES,CutOff=3 ){
    L.c <- list(
        KO= colnames(GES),Drug_Prism= colnames(DRS),Expr= colnames(Exp))
    CL.M <- Reduce(intersect,L.c)
    interactFeatures <- lapply(seq_len(nrow(Predtargets)), function(x){
    ## Drug response based on the drug ID
    DRS.f <- tryCatch(DRS[Predtargets[x,1],CL.M],error=function(e){NA})
    ## Gene effect scores based on the targeted gene
    GES.f <- tryCatch(GES[Predtargets[x,3],CL.M],error=function(e){NA})
    ## Groups based on the cut-off gene expression
    Exp.Group <- tryCatch(Exp[Predtargets[x,3],CL.M]< CutOff,
        error=function(e){NA})
    ## True is low, False is high.
    tryCatch(summary(lm(DRS.f ~ GES.f * Exp.Group))$coefficients[4,c(1,4)],
        error=function(e){NA})})
    names(interactFeatures) <- Predtargets$drugName
    interactFeatures
}


