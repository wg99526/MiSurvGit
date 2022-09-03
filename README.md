# MiCloudSGit

Title: MiCloud-S: An Integrative Web Cloud Platform for Microbiome Data Analysis with Survival Responses

Version: 1.0.0

Date: 2022-08-22

Maintainer: Won Gu <wpg5129@psu.edu> Hyojung Jang <hyojung.jang@stonybrook.edu> 

Description: MiCloud-S is an integrative web cloud platform for processing, analyzing and visualizing microbiome data with survival responses. MiCloud-S consists of a data processing module and following four data analytic modules as follows.

* **Data Processing** : Interactive procedures for 1) data inputs (.rdata, .rds, .biom, .txt, .csv, .tsv, .tre), 2) survival data and analytic plans (survival time, censored/event, follow-up period, subgroup analysis), 3) quality controls (kingdom, library size, mean proportion, taxonomic name), and 4) data transformations (alpha- and beta-diversity calculation, rarefaction, proportion, centered log-ratio, arcsine square root).

* **Module 1** : Comparative survival analysis between treatment groups, not involving microbiome data, with or without covariate adjustment.

* **Module 2** : Comparative analysis in microbial composition between treatment groups, not involving survival data, with or without covariate adjustment.

* **Module 3** : Association testing between microbial composition and survival responses with or without covariate adjustment.

* **Module 4** : Prediction modeling using microbial taxa at different taxonomic ranks on survival responses.


## URLs

* Web application (online implementation): https://223.194.200.99:3838
* GitHub repository (local implementation): https://github.com/wg99526/MiCloudSGit
 
## References

* Gu W, Koh H, Jang HJ, Lee B, Kang, B. MiCloud-S: An integrative web cloud platform for user-friendly microbiome data analysis with survival responses. (*_Submitted_*)

# Prerequites

shiny
```
install.packages("shiny")
```

# Launch App

```
library(shiny)

runGitHub("MiCloudSGit", "wg99526", ref = "main")
```

# Troubleshooting Tips

If you have any problems for using MiCloud-S, please report in Issues (https://github.com/wg99526/MiCloudSGit/issues) or email Won Gu (wpg5129@psu.edu) or Hyojung Jang (hyojung.jang@stonybrook.edu). 
