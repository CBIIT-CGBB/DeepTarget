\name{Depmap2DeepTarget}
\alias{Depmap2DeepTarget}
\title{
Retrieval and preparation of input data required from Depmap to Deeptarget package.
}
\description{
Retrieve gene expression, Cripr, mutation data from KO method, and drug matrix and then preperation the matrix compatible as input for Deeptarget.
}
\usage{
Depmap2DeepTarget(FileN,version)
}
\arguments{
    \item{FileN}{File Named used as input for DeepTarget: "CCLE_expression.csv", "CRISPRGeneEffect.csv"    ,"OmicsSomaticMutations.csv", or "secondary-screen-dose-response-curve-parameters.csv"}
    \item{version}{Version of data}
}
\value{
a data frame for each required input data
}
\author{
Trinh Nguyen, Ying Hu, and sanju
}

\examples{
library(readr)
library(depmap)
# expresion
CCLE.exp <- Depmap2DeepTarget("CCLE_expression.csv","19Q4")
}

