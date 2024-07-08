\name{OntargetM}
\alias{OntargetM}
\docType{data}
\title{An object containing a small part of the data from the Cancer Dependency Map (depmap.org) to demonstrate in DeepTarget pipeline}
\description{
An object containing Viability matrix after CRISPR-KO; Viability after Drug Treatment; Drug metadata from Broad, mutation matrix, and expression matrix with common cell-lines and common drugs. This is a subset of the total data due to memory constraints, full data can be downloaded from depmap.org/portal.
}
\usage{data("OntargetM")}
\format{A list of one dataframe and 4 matrices
    \describe{
        \item{\code{DrugMetadata}}{a dataframe containing 11 unique drugs as rownames with their associated information:  broad_id_trimmed as ID of the drug, name, target, drug_category, and moa as columns}
        \item{\code{secondary_prism}}{a viability scores matrix (after Drug Treatment) with 16 drugs as row names across 392 unique celllines as column names}
        \item{\code{avana_CRISPR}}{a Gene effect scores (after CRISPR-KO) matrix for 487 genes as row names across 392 unique celllines as column names}
        \item{\code{mutations_mat}}{Mutation binary matrix for 476 genes as row names across 392 unique cell lines as column names; 0 is WT; 1 is mutated}
        \item{\code{expression_20Q4}}{ Expression matrix for  550 genes as row names across 392 unique celllines as column names}
    }
}
\details{
For a full list data used in the paper, please use the link below to download data
}
\source{
DrugMetadata: Please download full data from this link \url{ https://depmap.org/repurposing/#:~:text=Corsello_supplemental_tables.xlsx}

Secondary prism: please download full data from this link \url{https://depmap.org/portal/download/all/?releasename=PRISM+Repurposing+19Q4&filename=secondary-screen-dose-response-curve-parameters.csv}

avana_CRISPR: please download full data from this link \url{https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q4&filename=CRISPRGeneEffect.csv}

mutations_mat:  Please download full data from this link \url{https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q4&filename=OmicsSomaticMutations.csv}

expression_20Q4: Please download full data of file named "CCLE_expression.csv" from this link \url{https://depmap.org/portal/download/all/}

}
\examples{
data(OntargetM)
}
\keyword{datasets}
