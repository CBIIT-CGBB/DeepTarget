
## perform the interaction based on cut-off

DoInteractExp <- function(Predtargets=Drug.sim,Exp=Exp,DRS=DRS, GES=GES,CutOff  = 3 )
  {
    L.c <- list( KO= colnames(GES),
                 Drug_Prism= colnames(DRS),
                 Expr= colnames(Exp))
    CL.M = Reduce(intersect,L.c)
  interactFeatures=lapply(1:nrow(Predtargets), function(x)
  {
    ## drug response based on the drug ID
    DRS.f=errHandle(DRS[Predtargets[x,1],CL.M])
    ## Gene effect scores based on the targeted gene
    GES.f=errHandle(GES[Predtargets[x,3],CL.M])
    ## groups based on the cut-off. Lower and higher based on the targed gene
    Exp.Group=errHandle(Exp[Predtargets[x,3],CL.M]< CutOff )
    ## True is low, False is high.

    errHandle(summary(lm(DRS.f ~
                             GES.f * Exp.Group))$coefficients[4,c(1,4)])

    }
)

  names(interactFeatures)=Predtargets$drugName
  interactFeatures
}


