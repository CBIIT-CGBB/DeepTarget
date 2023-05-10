## This is an example script for the package. It assumes that the data can contain multiple assays for the same drug.
rm(list=ls());
library(rstudioapi)
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path ))
library ( DeepTarget)
data ("OntargetM")
### Drugs used in the paper.
## can't locate the drug named GSK2830371 as the figure 5G.

drug.name <- c('atiprimod','AMG-232','pitavastatin','Ro-4987655','alexidine','RGFP966','dabrafenib','olaparib','CGM097','ibrutinib','palbociclib')
length(drug.name)
## Get the Drug ID
S.Drug <- OntargetM$DrugMetadata$broad_id_trimmed [which (OntargetM$DrugMetadata$name %in% drug.name)]
## Each of these drugs has two assays.
########
### this
length(S.Drug)
## 11
sec.prism.f <- OntargetM$secondary_prism[which ( row.names(OntargetM$secondary_prism) %in% S.Drug), ]
dim(sec.prism.f)
## 16
## Some of these drugs have multiple assays.

head(sec.prism.f)
KO.GES <- OntargetM$avana_CRISPR
dim(KO.GES)
dim(sec.prism.f)
## Calculate the similarity between these assays using the KO method.
List.sim <- NULL;
for ( i in 1:nrow(sec.prism.f)){
    DRS=as.data.frame(sec.prism.f[i,])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[i]
    out <- GetSim(row.names(sec.prism.f)[i],DRS=DRS, GES=KO.GES)
    List.sim [[length(List.sim) + 1]] <- out
}
names(List.sim) <- row.names(sec.prism.f)
dir.create( "./Result/")
saveRDS(List.sim,
        file = 'Result/similarity_KO_DrugTreatment.RDS')

###
length(List.sim)
####
metadata <- OntargetM$DrugMetadata
## Get the similarity for known gene targets from the drug (if there are multiple targeted genes, get the most similar gene)
DrugTargetSim <- PredTarget(Sim.GES.DRS=List.sim, D.M = metadata)
## Get the similarity for the gene which is most similar to the drug treatment.
DrugGeneMaxSim <- PredMaxSim(Sim.GES.DRS=List.sim, D.M = metadata)
## Mutant interaction and expression interaction with the known targeted gene.
d.mt <- OntargetM$mutations_mat
d.expr <- OntargetM$expression_20Q4
out.MutantTarget <- NULL;
out.LowexpTarget <- NULL;
### If duplicated exists:
for ( i in 1:nrow(sec.prism.f)){
    identical (row.names(sec.prism.f)  , DrugTargetSim[,1])
    DRS=as.data.frame(sec.prism.f[i,])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[i]
    ## For mutant:
    MutantInteract <- DoInteractMutant (Predtargets=DrugTargetSim[i,],Mutant=d.mt,DRS=DRS,GES=KO.GES)
    ## Assign the estimate as the strength and use the P value from the interaction model:
    TargetMutSpecificity <- data.frame(MaxTgt_Inter_Mut_strength=sapply(MutantInteract, function(x) x[1]), MaxTgt_Inter_Mut_Pval=sapply(MutantInteract, function(x) x[2]))
    out.MutantTarget <- rbind (out.MutantTarget,TargetMutSpecificity)
    ## For expression:
    ExpInteract <- DoInteractExp (DrugTargetSim[i,],d.expr,DRS=DRS,GES=KO.GES,CutOff = 2)
    ## Assign the estimate as the strength and use the P value from the interaction model.
    TargetExpSpecificity <- data.frame(MaxTgt_Inter_Exp_strength=sapply(ExpInteract, function(x) x[1]), MaxTgt_Inter_Exp_Pval=sapply(ExpInteract, function(x) x[2]))
    out.LowexpTarget <- rbind ( out.LowexpTarget,TargetExpSpecificity)
}

identical ( DrugTargetSim[,1],DrugGeneMaxSim[,1])

## Check whether the interaction is true or false based on a cut-off. Estimate and p value come from lm model
Whether_interaction_Ex_based= ifelse ( out.LowexpTarget$MaxTgt_Inter_Exp_strength <0 & out.LowexpTarget$MaxTgt_Inter_Exp_Pval <0.2,TRUE,FALSE)
## Mutation interaction with P <0.1
predicted_resistance_mutation = ifelse ( out.MutantTarget$MaxTgt_Inter_Mut_Pval<0.1,TRUE,FALSE)
### If desired, save how many cell lines have low expresion:
Pred.d <- cbind ( DrugTargetSim,DrugGeneMaxSim,out.MutantTarget,predicted_resistance_mutation, out.LowexpTarget,Whether_interaction_Ex_based)

Low.Exp = sapply(Pred.d[,3],function(x)errHandle(sum(d.expr[x,] < 2)) )
## Save this for later:
Pred.d$lowExpCount<-Low.Exp
##
### Obtain the similarity between the viablity scores from drug targeted treatment
## vs Gene effect score from KO method for the group that don't have primary targeted express.
## Prediction is in column 3.
## Only use the rows with no NA values.
idx <- which ( Pred.d$lowExpCount>0)
Pred.d.f <- Pred.d[idx ,]
Low.Exp.G = sapply(Pred.d.f[,3], function(x) errHandle(names(which(d.expr[x,]<2))))
identical ( names(Low.Exp.G),Pred.d.f[,3] )
## Only perform the test if the gene has at least some cell lines with low expression.

sim.LowExp <- NULL;
sec.prism.f.f <- sec.prism.f[idx,]
identical (row.names(sec.prism.f.f) ,Pred.d.f [,1])

for ( i in 1:nrow(Pred.d.f)){
    DRS.L= sec.prism.f.f[i,Low.Exp.G[[unlist(Pred.d.f[i,3])]]]
    DRS.L <- t(as.data.frame(DRS.L))
    row.names(DRS.L) <- Pred.d.f[i,1]
    out <- GetSim(Pred.d.f[i,1],DRS=DRS.L, GES=KO.GES)
    sim.LowExp [[length(sim.LowExp) + 1]] <- out
}
names(sim.LowExp) <-Pred.d.f[,1]
saveRDS(sim.LowExp,
        file = 'Result/similarity_KO_LowExp_DrugTreatment.RDS')

### Plot and show the top 5 genes that have the best correlation with the drug when the primary target is not expressed.

sim.LowExp.Strength=sapply(sim.LowExp, function(x) x[,2])
head(sim.LowExp.Strength)
sim.LowExp.Pval=sapply(sim.LowExp, function(x) x[,1])
head(sim.LowExp.Pval)
dim(sim.LowExp.Strength)
## We should create an rdata save it in the data folder if we can't create a good dummy file as the output.
## depending on whether we can add more data in the pacakge.
## the best is to save both similary for all and for low.
pdf("Result/sim.low.exp.plot.pdf")
par(mar=c(4,4,5,2), xpd=TRUE, mfrow=c(2,2));
plotSim (dx=sim.LowExp.Pval,dy=sim.LowExp.Strength,clr=colorRampPalette(c("lightblue",'darkblue')), plot=TRUE)
dev.off();

## Record the top 5 genes to the pred object.

L.topG <- NULL;
for ( i in 1:ncol(sim.LowExp.Strength)){
    top.5 <- names (sort(sim.LowExp.Strength[,i], decreasing=TRUE)[1:5])
    top.5  <- unlist(top.5);
    top.5    <- paste(top.5, collapse=" ");
    L.topG <- rbind ( L.topG, top.5 )
}
Pred.d$top5GeneWlowEx <- "NA"
Pred.d$top5GeneWlowEx [idx] <- L.topG

H.topG <- NULL;
simExp.Strength <- sapply(List.sim, function(x) x[,2])
dim(simExp.Strength)
for ( i in 1:ncol(simExp.Strength)){
    top.5 <- names (sort(simExp.Strength[,i], decreasing=TRUE)[1:5])
    top.5  <- unlist(top.5);
    top.5    <- paste(top.5, collapse=" ");
    H.topG <- rbind ( H.topG, top.5 )
}

Pred.d$top5Gene_Ex <- H.topG

write.csv (Pred.d,"./Result/Prediction_sim_KO_DrugTreatment.csv" )
## let plot for mutation =TRUE

## plot
DOI = 'dabrafenib'
GOI ='BRAF'
head(Pred.d)
## The first two columns are the drugs having two assays.
## We know that sec.prism.f has the same order as Pred.d.
which.mut <- which (Pred.d$predicted_resistance_mutation==TRUE);
## Plot mutants:
### Preparing data:
for ( i in 1:length(which.mut)){
    cr.i <- which.mut[i]
    DOI = Pred.d[cr.i,2]
    GOI =Pred.d[cr.i,3]
    DRS=as.data.frame(sec.prism.f[cr.i,])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[cr.i]
    head(DRS)
    pdf ( paste0("./Result/mutant_", row.names(Pred.d)[cr.i],".pdf"))
    out <- DMB (DN=DOI,GN=GOI,Pred=Pred.d[cr.i,],Mutant=d.mt,DRS= DRS,GES= KO.GES,plot=TRUE)
    print (out)
    dev.off();
}

## Primary targets:
which.exp <- which (Pred.d$Whether_interaction_Ex_based==FALSE);
for ( i in 1:length(which.exp )){
    cr.i <- which.exp[i]
    DOI = Pred.d[cr.i,2]
    GOI =Pred.d[cr.i,'MaxTargetName']
    DRS=as.data.frame(sec.prism.f[cr.i,])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[cr.i]
    head(DRS)
    pdf ( paste0("./Result/expr_", GOI,row.names(Pred.d)[cr.i],".pdf"))
    out <- DTR ( DN=DOI,GN=GOI,Pred=Pred.d[cr.i,], Exp=d.expr,DRS= DRS,GES=KO.GES,CutOff= 2)

    print (out)
    dev.off();
}

## Based on the secondary targets:
which.exp <- which (Pred.d$Whether_interaction_Ex_based==FALSE);
dim(Pred.d)
for ( i in 1:length(which.exp )){
    cr.i <- which.exp[i]
    DOI = Pred.d[cr.i,2]
    GOI =Pred.d[cr.i,'BestTargetGene']
    DRS=as.data.frame(sec.prism.f[cr.i,])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[cr.i]
    head(DRS)
    pdf ( paste0("./Result/expr_", GOI,row.names(Pred.d)[cr.i],".pdf"))
    out <- DTR ( DN=DOI,GN=GOI,Pred=Pred.d[cr.i,], Exp=d.expr,DRS= DRS,GES=KO.GES,CutOff= 2)
    print (out)
    dev.off();
}
## Plot the correlation for the predited target.
DOI = 'atiprimod'
GOI = 'SOX10'
#GOI = 'MITF'
idx <- which( Pred.d[,2]==DOI)
dim(sec.prism.f)
identical ( row.names(sec.prism.f) , Pred.d$DrugID)
for ( i in 1:length(idx)){
    DRS=as.data.frame(sec.prism.f[idx[i],])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[idx[i]]
    head(DRS)
    pdf ( paste0("./Result/Cor_plot_Predicted_Of", GOI,'of_',row.names(Pred.d[idx[i],]),".pdf"));
    plotCor(DN=DOI,GN=GOI,Pred=Pred.d[idx[i],],DRS= DRS,GES= KO.GES,plot=TRUE);
    dev.off()
}



