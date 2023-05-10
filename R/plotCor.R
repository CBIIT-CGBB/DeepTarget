
plotCor <- function(DN=Drugname,GN=GOI,Pred=Pred,DRS= DRS,GES= GES,plot=TRUE){
    ## need to add to catch error.
    ## Checking to see whether the drug name exists in the Pred object.
    ## Get the overlapping drugs

    if (nrow (Pred)==0) {stop("Pred should contain drug of interest")}
    L.c <- list( KO= colnames(GES),
                 Drug_Prism= colnames(DRS))
    CL.M = Reduce(intersect,L.c)
    if ( DN %in% Pred[,2]){
        if (plot) {
            dat.c <- data.frame(GES.c = GES [GN,CL.M], DRS.c = DRS[Pred[DN,1],CL.M])
            ggplot(dat.c, aes(x = GES.c, y = DRS.c))+
                #scale_color_manual(values = c("cyan",'red')) +
                geom_point()+
                stat_smooth(method = 'lm')+
                theme_bw(base_size = 15)+
                labs(x = paste0(GN,'_Viability after CRISPR-KO'), y = paste0('Response to ', Pred[DN,2]))+
                stat_cor(label.y = c(1.3,1.4))
        }
    }else{
        stop("The drug of interest does not exist. Please double-check with the target prediction output.")
    }
}
