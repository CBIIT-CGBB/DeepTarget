# DeepTarget
## About
DeepTarget is a computational tool for deep characterization of cancer drugs’ MOA by integrating existing large-scale genetic and drug screens. Spanning ~1500 drugs across ~18K possible target genes, DeepTarget predicts: (1) a drug’s primary target(s), (2) whether it specifically targets the wild-type or mutated target forms, and (3) the secondary target(s) that mediate its response when the primary target is not expressed. We first tested and successfully validated DeepTarget in a total of eleven unseen gold-standard datasets, with an average AUC of 0.82, 0.77, and 0.92 for the above three prediction abilities, respectively. We then proceed to use it in a wide range of applications: First, we find that DeepTarget’s predicted specificity of a drug to its target is strongly associated with its success in clinical trials and is higher in its FDA-approved indications. Second, DeepTarget predicts candidate drugs for targeting currently undruggable cancer oncogenes and their mutant forms. Finally, DeepTarget predicts new targets for drugs with unknown MOA and new secondary targets of approved drugs. Taken together,DeepTarget is a new computational framework for accelerating drug prioritization and its target discovery by leveraging large-scale genetic and drug screens.

Full artice (preprint): https://www.biorxiv.org/content/biorxiv/early/2022/10/19/2022.10.17.512424.full.pdf

# Quick Installation
## Install devtools in your local R environment.
install.packages("devtools")
## Call devtools
library(devtools)
## Install DeepTarget from the GitHub
devtools::install_github("CBIIT-CGBB/DeepTarget")


Please refer to the R script in Example/ExWorkflow.R for a demonstration of the code.
