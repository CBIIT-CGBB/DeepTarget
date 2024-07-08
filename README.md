
# DeepTarget
`DeepTarget` performs deep characterization of cancer drugs mechanism of action by integrating large-scale genetic and drug screens from the Cancer Dependency Map project run by the Broad Institute (depmap.org)  and aims to find drug's primary target(s), predict whether the drug specifically targets the wild-type or mutated target forms, and predict the secondary target(s) that mediate its response when the primary target is not expressed.

## Installation
To install depmap, the BiocManager Bioconductor Project Package Manager is required. If BiocManager is not already installed, it will need to be done so beforehand with (within R) install.packages("BiocManager"). Once it is installed, DeepTarget can be installed from Biocondctor:
``` r
install.packages("BiocManager")
BiocManager::install("DeepTarget")
```
### To install the version from GitHub, use

``` r
install.packages("BiocManager")
BiocManager::install("CBIIT-CGBB/DeepTarget")
```

## Example

Please refer to the file named DeepTarget_Vignette.Rmd in the directory vignettes for a demonstration of how the package can be used.

## Data
For the purpose of demonstration, we include an OntargetM object that contains a subset of data from the Cancer Dependency Map (depmap.org). The full data set it too large for package memory constraints. For details and links to download the full data set, please use the  `??OntargetM` command.

## Citation
If you use `DeepTarget`, please consider adding the following
citation from bioRxiv preprint: (https://www.biorxiv.org/content/10.1101/2022.10.17.512424v1)
@article {Sinha2022.10.17.512424,
	author = {Sanju Sinha and Neelam Sinha and Eytan Ruppin},
	title = {Deep characterization of cancer drugs mechanism of action by integrating large-scale genetic and drug screens},
	elocation-id = {2022.10.17.512424},
	year = {2022},
	doi = {10.1101/2022.10.17.512424},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2022/10/19/2022.10.17.512424},
	eprint = {https://www.biorxiv.org/content/early/2022/10/19/2022.10.17.512424.full.pdf},
	journal = {bioRxiv}
}
