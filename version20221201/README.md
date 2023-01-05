## DeepTarget: Predicting in-depth Mechanism of Cancer Drugs
 <img width="605" alt="Screen Shot 2022-12-01 at 11 51 50 AM" src="https://user-images.githubusercontent.com/29957280/205112105-64f2c923-58c5-4d99-8cf6-09b88ee9e0e7.png">


### Prerequisites
> Software: RStudio, R

> Packages required: parallel,tictoc,fgsea, ggplot2, ggpubr, statar, interactions,stats, reshape2, pROC, data.table,ggrepel, grid,cowplot,gridExtra, stringb, DEGreport, readxl, stringr

> Please download the data required to execute the scripts from here: https://doi.org/10.5281/zenodo.7367726
Following this, please move them to a folder name 'Data'. Please move 'Data' to the repository 'drug.OnTarget' & set that as your 'working_dir' in the scripts.

### Core Steps of DeepTarget:
> Tools >> 'Core Steps'

> Step0_Write_Functions.Rmd
Define all function required for the analysis and validation of our pipeline. 

> Step0B_Data_and_Libraries.Rmd
Pre-loading and preprocessing data and libraries required.
* Required File:
 * onTarget (onTarget_v3.0.RDS): A single object is a curation of all the DepMap data including: viability screens after CRISPR-Cas9 & drug treatment, expression & mutation of cell lines.
 * KnowTraget_prediction (KnownTarget_predictions_v4.RDS): Metadata of Drugs and their target+DeepTarget's top predicted target information 
 * drugVScrispr_corr_features_list (drugVScrispr_corr_features_list.RDS): DKS Score & P for each gene vs each drug. 
 * drugVScrispr_corr_features_list_secondary (drugVScrispr_corr_features_list_secondary.RDS): Secondary DKS Score & P for each gene vs each drug.
 * human_Kinases (human_Kinases.xlsx): List of all Human Kinases genes

---
### The result of the below scripts are saved in the Data folder.

1) Step 1: Step1_Create_Drug_KO_similarity_matrix.Rmd:
> Step_0A_data_curation_and saving.Rmd
* Input:
  * onTarget (onTarget_v3.0.RDS)
* Output:
  * drugVScrispr_corr_features_list (drugVScrispr_corr_features_list.0.RDS)

In this Step, we compared large-scale drug response screens (PRISM) with ~1500 cancer drugs & genome-wide CRISPR-Cas9 knockout viability profiles (CRISPR-KO essentiality, AVANA) from DepMap, which were commonly performed across 371 cancer cell lines. We calculated a (drug, target) pair-wise similarity score (Pearson Correlation) between viability after drug treatment and that observed for each gene CRISPR-KO termed it the Drug-KO Similarity score (DKS score & P; Object Name: drugVScrispr_corr_features_list). This score can range from -1 to 1. Using this DKS Score, we curated data KnownTarget_predition which will be majorly used the futher analysis and validation of our Pipeline. The result of the correlation is saved in the Data folder.

---
2) Step 2: Step2_Retreive_DKS_Score_KnownTargets.Rmd
* Input:
  * onTarget (onTarget_v3.0.RDS)
* Output:
  * DKS Score & P for known Target of each drug & best Predicted Target by DeepTarget (KnownTarget_predictions.RDS)

Retreive a DKS Score for Known Target of a drug; also provide a best predicted Primary Target for each drug.
---

3) Step 3: Step3_Primary_Target_at_Pathway_Resolution.Rmd
* Input:
  * onTarget (onTarget_v3.0.RDS) 
* Output:
  * Secondary DKS Score & P for each pathway vs each drug.

Provide each pathway a likelihood score of being a secondary Mechanism of Action.

---
4) Step 4: Step3_Find_Secondary_Targets.Rmd.
* Input:
  * onTarget (onTarget_v3.0.RDS) 
* Output:
  * Secondary DKS Score & P for each gene vs each drug (drugVScrispr_corr_features_list_secondary.RDS)

Provide each gene a likelihood score of being a secondary target of drug. (Exploratory: We first hypothesize and tested whether the Drugs where the primary target is not expressed in most cell lines, whether the Primary Target score of Known Target very low. Second, We test drugs whose target KO and drug response viability similarity dimnishes when target expression is not expressed.) Crux: In this step, we repeated our target identification process (Step 1) in cell lines where their primary targets are not expressed to produce a secondary-DKS Score for such drugs. This saved dataset is used in secondary target analysis and figure generation.

---

5) Step 5: step4_Mutant_Specificity.Rmd.
* Input:
  * onTarget (onTarget_v3.0.RDS) 
* Output:
  * Mutant Specificity Score & P for each drug vs their known target (Target_Mutation_specificity.RDS)

  Here, we compute the Mutant Specificity Score (MS Score) of a drug to its known target: We reason that if a drug more specifically targets a mutant form of protein, in the cell lines with this mutant form, the similarity between viability after drug treatment and target CRISPR-KO (their DKS score) would be significantly higher than in the cell lines with WT protein. We mathematically model this dependency of the DKS score on the mutation status of the target by employing a regression between the drug and the target’s response across the cell lines, including the mutation status as an interaction term. We call this score the mutant-specificity score (this composes the second output of DeepTarget), where a positive mutant-specificity score indicates that a drug differentially targets the mutant form of the target protein more than WT and vice-versa for the negative score.
---
