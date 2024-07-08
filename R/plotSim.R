## Function to plot the similarity based on correlation values and p values.
plotSim <- function(dx,dy,clr=NULL, plot=TRUE){
    if (is.null(clr)){ clr <- colorRampPalette(c("lightblue",'darkblue'))}
    dx <- as.matrix(dx)
    dy <- as.matrix(dy)
    if (plot) {
        if (identical(dimnames(dx), dimnames(dy))){
            for ( i in seq_len(ncol(dx))){
                x <- -log10 (dx[,i])
                y <- dy[,i]
                top.5 <- names (sort(dy[,i], decreasing=TRUE)[seq_len(5)])
                y.clr <- clr(10)[as.numeric(cut(y,breaks = 10))]
                x1   <- (min(x) + max(x))/2 - 2;
                x2   <- (min(x) + max(x))/2 + 2;
                y1   <- max(y) + (min(abs(y)) + max(abs(y))) * 0.40;
                y2   <- max(y) + (min(abs(y)) + max(abs(y))) * 0.55;
                plot(x, y, pch = 19, col = y.clr , 
                    xlab="-log10(Pval)",
                    ylab= paste0("correlation strength ",colnames(dx)[i]));
                text(x[top.5],  y [top.5],  top.5,cex=0.65, pos=2.8,col="black")
                colorBar(x1, y1, x2, y2, v.min=min(y),
                    v.max=max(y),col.l=clr(100))
            }
        }else{
            stop("the dimension of input is not the same, please double check")}
    }
}
