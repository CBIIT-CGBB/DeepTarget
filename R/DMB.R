## Function to calculate Drug Mutant Binding (DMB)
DMB <- function(DrugName,GOI,Pred,Mutant,DRS,GES,plot=TRUE){
    if (nrow (Pred)==0) {stop("Pred should contain drug of interest")}
    L.c <- list(
        KO= colnames(GES),Drug_Prism= colnames(DRS),Mutation= colnames(Mutant))
    CL.M <- Reduce(intersect,L.c)
    if ( DrugName %in% Pred[,2]){
        if (plot) {
            mutant  <- factor(Mutant[GOI,CL.M]);
            colcode <- mutant
            ## zero WT, 1 mutant.
            levels(colcode ) <- c('WT', 'Mutant')
            GES.c <- GES[GOI,CL.M]
            DRS.c <- DRS[Pred[DrugName,1],CL.M]
            Mutant.c <- colcode
            dat.c <- data.frame(GES.c,DRS.c,Mutant.c)
            ggplot(dat.c,
                aes(x = GES.c, y = DRS.c, color = factor(dat.c$Mutant.c)))+
                stat_smooth(method = 'lm')+
                theme_bw(base_size = 15)+
                labs(
                    x = 'Viability after CRISPR-KO', 
                    y = paste0('Response to ',Pred[DrugName,2]),
                    color =  paste0( GOI,'\nMutation\nstatus'))+
                stat_cor(label.y = c(1.3,1.4))
            }
    }else{
        stop("The drug of interest is not existed. Please double check.")}
}


