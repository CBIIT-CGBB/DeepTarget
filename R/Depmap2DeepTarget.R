## This version is pull the data from publicly available files that were shared by the depmap consortium on figshare. 
Depmap2DeepTarget <- function(FileN,version){
    out <- depmap::dmsets() |> dplyr::filter(grepl(version, title))
    if (nrow(out)==0){stop("please check FileName and its version")}
    for ( i in seq_len(nrow(out))){
        d.id <- out$dataset_id[i]
        Found.i <-  depmap::dmfiles() |>filter(dataset_id == d.id) |> filter(name == FileN)
        if (nrow(Found.i)>0){break}
    }
    if (nrow(Found.i)==0){
        stop("This file is not available for this version")
    }else{
        ## download and read the file.
        out.f <- depmap::dmfiles() |>filter(dataset_id == d.id) |> filter(name == FileN) |>
            dmget() |>
            read_csv()
        if (FileN=="CCLE_expression.csv" || FileN=="CRISPRGeneEffect.csv"){
            out.t <- t(out.f)
            colnames(out.t) <-  out.t[1,]
            out.t <- out.t[-1,]
            out.d <- as.data.frame(out.t)
            out.d$Gene <- row.names(out.d)
            # pre-processed
            Gene_name <- vapply(strsplit(out.d$Gene, '[ ()]'), "[", 1,FUN.VALUE = as.character(length(out.d$Gene)) )
            out.m <- cbind (Gene_name,out.d)
        }
        if (FileN=="OmicsSomaticMutations.csv"){
            mu_frame <- data.frame(ColumnName=out.f$DepMap_ID,
            RowName=out.f$HugoSymbol)
            ## remove silent one.
            filter_var_class <- "SILENT"
            f.i <- which(out.f$VariantInfo %in% filter_var_class)
            mu_frame <- mu_frame[-f.i,]
            # Unique column and row names
            unique_cols <- unique(mu_frame$ColumnName)
            unique_rows <- unique(mu_frame$RowName)
            # Create an empty matrix
            out.m <- matrix(0, nrow = length(unique_rows), ncol = length(unique_cols),
                dimnames = list(unique_rows, unique_cols))
            # Fill the matrix
            for(j in seq_len(nrow(out.m))){
                out.m[mu_frame$RowName[j], mu_frame$ColumnName[j]] <- 1
            }
            out.m <- as.data.frame(out.m)}
        if (FileN=="secondary-screen-dose-response-curve-parameters.csv"){
            # Area under the curve(AUC) is used; also, retrieve meta data for later use.
            out.f$Broad_id_trimmed <- vapply(strsplit(out.f$broad_id, "[-]"), "[", 2, FUN.VALUE = as.character(length(out.f$broad_id)))
            out.d.Auc <- data.frame(ColumnName=out.f$depmap_id,RowName=out.f$Broad_id_trimmed,Value=out.f$auc)
            out.d.Auc.m <- xtabs(Value ~ RowName + ColumnName, data = out.d.Auc)
            attr(out.d.Auc.m, "class") <- NULL
            attr(out.d.Auc.m, "call") <- NULL
            # annotation
            out.f.anno <- out.f[,c("broad_id","Broad_id_trimmed","name","target", "moa","smiles")]
            out.d.anno <- as.data.frame(out.f.anno)
            out.d.anno <- out.d.anno[!duplicated(out.d.anno$Broad_id_trimmed),]
            row.names(out.d.anno) <- out.d.anno$Broad_id_trimmed
            out.d.anno <- out.d.anno[row.names(out.d.Auc.m),]
            out.m <- cbind ( out.d.anno,out.d.Auc.m)}
            # return the dataframe if found
        return(out.m)
    }
}


