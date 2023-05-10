
DMB <- function(DN=Drugname,GN=GOI,Pred=Pred,Mutant=Mutant,DRS= DRS,GES= GES,plot=TRUE){
    ## need to add to catch error.
    ## checking to see whether the drug name exist in the Pred object.
    ## get the overlapped.

    if (nrow (Pred)==0) {stop("Pred should contain drug of interest")}
    L.c <- list( KO= colnames(GES),
                 Drug_Prism= colnames(DRS),
                 Mutation= colnames(Mutant))
    CL.M = Reduce(intersect,L.c)
    if ( DN %in% Pred[,2]){
        if (plot) {
            mutant  <- factor(Mutant[GN,CL.M]);
            colcode <- mutant
            ## zero WT, 1 mutant.
            levels(colcode ) <- c('WT', 'Mutant')
            dat.c <- data.frame(GES.c = GES [GN,CL.M], DRS.c = DRS[Pred[DN,1],CL.M], Mutant.c = colcode)
            ggplot(dat.c, aes(x = GES.c, y = DRS.c, color = factor(dat.c$Mutant.c)))+
                #scale_color_manual(values = c("cyan",'red')) +
                # geom_point()+
                stat_smooth(method = 'lm')+
                theme_bw(base_size = 15)+
                labs(x = 'Viability after CRISPR-KO', y = paste0('Response to ', Pred[DN,2]), color =  paste0( GN,'\nMutation\nstatus'))+
                stat_cor(label.y = c(1.3,1.4))
        }
    }else{
        stop("The drug of interest does not exist. Please double-check with target prediction output.")
    }
}


