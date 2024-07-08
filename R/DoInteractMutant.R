## Function to calculate the interaction for a mutant target
DoInteractMutant <- function(Predtargets,Mutant,DRS,GES){
    # Retrieve only the overlapping samples
    L.c <- list( KO= colnames(GES),
                Drug_Prism= colnames(DRS),
                Mutation= colnames(Mutant))
    CL.M <- Reduce(intersect,L.c)
    interactFeatures <- lapply(seq_len(nrow(Predtargets)), function(x){
            ##  Drug response scores (DRS) based on the drug ID
            DRS.f <- tryCatch(DRS[Predtargets[x,1],CL.M],error=function(e){NA})
            ## Gene effect scores (GES) based on the targeted gene
            GES.f <- tryCatch(GES[Predtargets[x,3],CL.M],
                error=function(e){NA})
            Mutant.Group <- tryCatch(Mutant[Predtargets[x,3],CL.M] )
            ## True is low, False is high.
            tryCatch(summary(
                lm(DRS.f ~ GES.f * Mutant.Group))$coefficients[4,c(1,4)],
                error=function(e){NA})
            })
    names(interactFeatures) <- Predtargets$drugName
    interactFeatures
}


