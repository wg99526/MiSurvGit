# MiSurv

Title: MiSurv: An Integrative Web Cloud Platform for Comprehensive Microbiome Data Analysis with Survival Responses

Version: 2.0.1

Date: 2023-07-12

Maintainer: Won Gu <wpg5129@psu.edu> Hyojung Jang <hyojung.jang@northwestern.edu> 

Description: MiSurv is an integrative web cloud platform for processing, analyzing and visualizing microbiome data with survival responses. MiSurv consists of a data processing module and its following four data analytic modules as below.

* **Data Processing** : Interactive procedures for (1) data inputs (.rdata, .rds, .biom, .txt, .csv, .tsv, .tre), (2) survival data and analytic plans (survival time, censored/event, follow-up period, subgroup analysis), (3) quality controls (kingdom, library size, mean proportion, taxonomic name), and (4) data transformations (alpha- and beta-diversity calculation, rarefaction, proportion, centered log-ratio, arcsine square root).

* **Module 1** : Comparative survival analysis between treatment groups, not involving microbiome data, with or without covariate adjustment(s).

* **Module 2** : Comparative analysis in microbial composition between treatment groups, not involving survival data, with or without covariate adjustment(s).

* **Module 3** : Association testing between microbial composition and survival responses with or without covariate adjustment(s).

* **Module 4** : Prediction modeling using microbial taxa at different taxonomic ranks on survival responses.


## URLs

* Web application (online implementation): http://misurv.micloud.kr
* GitHub repository (local implementation): https://github.com/wg99526/MiSurvGit
 
## Example Data

* We stored the all the final processed data that we used in our paper (Gu et al., 2023) in the above Data folder as ‘phyloseq’ format (see the file named ‘biom.Rdata’) and also as four individual files (see the files named ‘otu.tab.txt’, ‘sam.dat.txt’, ‘tax.tab.txt’, and ‘tree.tre’).
* The raw sequence data are deposited in QIITA (https://qiita.ucsd.edu) with the ID number 10508 (https://qiita.ucsd.edu/study/description/10508). 

## References

* Gu, W., Koh, H., Jang, H., Lee, B., Kang, B. (2023). MiSurv: An integrative web cloud platform for user-friendly microbiome data analysis with survival responses. *_Microbiology Spectrum_* 11(3): e05059-22

# Prerequites

* For local implementation using this GitHub repository, recent version of **R (>4.2.0)** is needed.

shiny
```
install.packages("shiny")
```

# Launch App

```
library(shiny)

runGitHub("MiSurvGit", "wg99526", ref = "main")
```

# Troubleshooting Tips

If you have any problems for using MiSurv, please report in Issues (https://github.com/wg99526/MiSurvGit/issues) or email Won Gu (wpg5129@psu.edu) or Hyojung Jang (hyojung.jang@northwestern.edu).
