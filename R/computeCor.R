# Function to compute the correlation for a drug.
computeCor <- function(DrugName,DRS,GES){
    common.c <- intersect(colnames(DRS),colnames(GES))
    DRS.c<- DRS[DrugName,common.c ]
    GES.c<- GES[,common.c]
    out.Sim.GES.DRS <- vapply(
        seq_len(nrow(GES.c)),
        function(x) unlist(cor.test(DRS.c, GES.c[x,])),FUN.VALUE=character(10))
    p.value <- t(out.Sim.GES.DRS["p.value",])
    estimate.cor <- t(out.Sim.GES.DRS["estimate.cor",])
    FDR <- p.adjust(p.value,method = 'fdr')
    out.Sim.GES.DRS.c <- cbind(as.numeric(p.value),as.numeric(estimate.cor),as.numeric(FDR))
    row.names(out.Sim.GES.DRS.c) <- row.names(GES.c)
    colnames(out.Sim.GES.DRS.c) <- c("p.value","estimate.cor","FDR")
    out.Sim.GES.DRS.c
}
