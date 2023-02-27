rm(list = ls())

list.of.packages <- c('seqinr', 'shinydashboard', 'dashboardthemes', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable', 
                      'DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'zCompositions', 'dplyr', 'forestplot', 'quantreg', 'fossil', 'picante',
                      'entropart', 'lme4', 'lmerTest', 'broom.mixed', 'gee', 'geepack', 'dirmult', 'robustbase', 'robCompositions', 'BiasedUrn',
                      'CompQuadForm', 'GUniFrac', 'ecodist', 'MiRKAT', 'gridExtra', 'ggplot2', 'patchwork', 'ggthemes', 'erer', 'DiagrammeR', 'stringr',
                      'devtools', 'betareg', 'nlme', 'glmm', 'remotes', 'gridGraphics', 'compositions', 'rgl', 'vegan3d', 'pca3d', 'jpeg', 'splitTools', 
                      'survival', 'survminer', 'coin', 'randomForestSRC', 'kableExtra', 'caret', 'randomForest', 'glmnet')

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!require('phyloseq')) remotes::install_github('joey711/phyloseq')
if(!require('biomformat')) remotes::install_github('joey711/biomformat')


library(seqinr)
library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(tidyverse)
library(phyloseq)
library(plotly)
library(shinyWidgets)
library(shinyjs)
library(googleVis)
library(xtable)
library(DT)
library(htmltools)
library(biomformat)
library(phangorn)
library(bios2mds)
library(zip)
library(randomForestSRC)
library(rgl)
library(vegan3d)
library(pca3d)
library(jpeg)
library(splitTools)

source("Source/MiDataProc.Data.Upload.R")
source("Source/MiDataProc.Alpha.Cross.Sectional.R")
source("Source/MiDataProc.Beta.Cross.Sectional.R")
source("Source/MiDataProc.Taxa.Cross.Sectional.R")
source("Source/MiDataProc.Surv.Model1.R")
source("Source/MiDataProc.Surv.Model3.Alpha.R")
source("Source/MiDataProc.Surv.Model3.Beta.R")
source("Source/MiDataProc.Surv.Model3.Taxa.R")
source("Source/MiDataProc.Surv.Model4.R")

# COMMENTS
{
  TITLE = p("MiSurv: An Integrative Web Cloud Platform for Microbiome Data Analysis with Survival Responses", style = "font-size:18pt")
  HOME_COMMENT = p(strong("MiSurv", style = "font-size:15pt"), "is an integrative web cloud platform for processing, analyzing and visualizing microbiome data with survival responses. MiSurv consists of a data processing module and its following four data analytic modules: (1) Module 1: Comparative survival analysis between treatment groups, (2) Module 2: Comparative analysis in microbial composition between treatment groups, (3) Module 3: Association testing between microbial composition and survival responses, (4) Module 4: Prediction modeling using microbial taxa on survival responses. More details are as follows.", style = "font-size:13pt")
  
  HOME_COMMENT1 = p(strong("Data Processing : "), "Interactive procedures for (1) data inputs (.rdata, .rds, .biom, .txt, .csv, .tsv, .tre), (2) survival data and analytic plans (survival time, censored/event, follow-up period, subgroup analysis), (3) quality controls (kingdom, library size, mean proportion, taxonomic name), and (4) data transformations (alpha- and beta-diversity calculation, rarefaction, proportion, centered log-ratio, arcsine square root).", style = "font-size:13pt")
  HOME_COMMENT2 = p(strong("Module 1:"), "Comparative survival analysis between treatment groups, not involving microbiome data, with or without covariate adjustment(s). ", style = "font-size:13pt")
  HOME_COMMENT4 = p(strong("Module 3:"), "Association testing between microbial composition and survival responses with or without covariate adjustment(s).", style = "font-size:13pt")
  HOME_COMMENT3 = p(strong("Module 2:"), "Comparative analysis in microbial composition between treatment groups, not involving survival data, with or without covariate adjustment(s).", style = "font-size:13pt")
  HOME_COMMENT5 = p(strong("Module 4:"), "Prediction modeling using microbial taxa at different taxonomic ranks on survival responses", style = "font-size:13pt")
  
  HOME_COMMENT6 = p(strong("URLs:"), " Web server (online implementation):", tags$a(href = "http://misurv.micloud.kr", "http://misurv.micloud.kr"), 
                    "; GitHub repository (local implementation):", 
                    tags$a(href = "https://github.com/wg99526/misurvgit", "https://github.com/wg99526/misurvgit"), style = "font-size:13pt")
  HOME_COMMENT7 = p(strong("Maintainers:"), " Won Gu (", tags$a(href = "wpg5129@psu.edu", "wpg5129@psu.edu"), 
                    "); Hyojung Jang (", tags$a(href = "hyojung.jang@stonybrook.edu", "hyojung.jang@stonybrook.edu"), ")", style = "font-size:13pt")
  HOME_COMMENT8 = p(strong("Reference:"), "Gu W, Koh H, Jang HJ, Lee B, Kang B. MiSurv: An integrative web cloud platform for user-friendly microbiome data analysis with survival responses. (in review)", style = "font-size:13pt")
  
  
  
  INPUT_PHYLOSEQ_COMMENT1 = p("Description:", br(), br(), "This should be an '.Rdata' or '.rds' file, and the data should be in the 'phyloseq' format (see", tags$a(href = "https://bioconductor.org/packages/release/bioc/html/phyloseq.html", "https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), ") The phyloseq object should contain all the four necessary data, feature (OTU or ASV) table, taxonomic table, metadata/sample information, and phylogenetic tree.", 
                              br(), br(), "Details:", br(), br(), 
                              "1) The feature table should contain counts, where rows are features (OTUs or ASVs) and columns are subjects 
                              (row names are feature IDs and column names are subject IDs).", br(),
                              "2) The taxonomic table should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                              (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                              'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(),
                              "3) The metadata/sample information should contain variables for the subjects about host phenotypes, medical interventions, 
                              disease status or environmental/behavioral factors, where rows are subjects and columns are variables 
                              (row names are subject IDs, and column names are variable names).", br(), 
                              "4) The phylogenetic tree should be a rooted tree. Otherwise, MiSurv automatically roots the tree through midpoint rooting (phangorn::midpoint). 
                              The tip labels of the phylogenetic tree are feature IDs.", br(), br(), 
                              "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between feature table and metadata/sample information. 
                              MiSurv will analyze only the matched features and subjects."
                              , style = "font-size:11pt")
  
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download example microbiome data 'biom.Rdata' in the 'phyloseq' format.", br(), 
  "For more details about 'phyloseq', see", 
                              tags$a(href = "https://bioconductor.org/packages/release/bioc/html/phyloseq.html", "https://bioconductor.org/packages/release/bioc/html/phyloseq.html"),
                              br(),
  p(strong("Data description:"), "This example data are the public gut microbiome data (Zhang et al. 2018) we used in our real data applications (Gu et al., in review). The raw sequence data are deposited in QIITA (", tags$a(href = "https://qiita.ucsd.edu", "https://qiita.ucsd.edu"), ") with the ID number 10508 (", tags$a(href = "https://qiita.ucsd.edu/study/description/10508", "https://qiita.ucsd.edu/study/description/10508"), "). More detailed sample extraction and raw sequence data processing procedures are described in (Zhang et al. 2018; Gu et al., in review). "), 
  
  br(), 
                              "> setwd('/yourdatadirectory/')", br(), br(), 
                              "> load(file = 'biom.Rdata')", br(), br(), 
                              "> library(phyloseq)", br(), br(), 
                              " > otu.tab <- otu_table(biom)", br(), 
                              " > tax.tab <- tax_table(biom)", br(), 
                              " > tree <- phy_tree(biom)", br(), 
                              " > sam.dat <- sample_data(biom)", br(), br(), 
                              " > sam.dat$T1Dweek",   HTML('&emsp;'),HTML('&nbsp;'), HTML('&nbsp;'),  "# Survival Time", br(), 
                              " > sam.dat$T1D ", HTML('&emsp;'),HTML('&emsp;'),HTML('&emsp;'), HTML('&nbsp;'), HTML('&nbsp;'), "# Event", br(), 
                              " > sam.dat$Antibiotics ", HTML('&emsp;'), HTML('&nbsp;'), "# Treatment", br(), br(),
                              "You can check if the features are matched and identical across feature table, taxonomic table and 
                              phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information 
                              using following code.", br(), br(), 
                              " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                              " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT = p("Description:", br(), br(), 
                                    "1) The feature table (.txt or .csv) should contain counts, where rows are features (OTUs or ASVs) and columns are subjects 
                                    (row names are feature IDs and column names are subject IDs). Alternatively, you can upload .biom file processed by QIIME.", br(), 
                                    "2) The taxonomic table (.txt) should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                                    (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                                    'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'). Alternatively, you can upload .tsv file processed by QIIME.", br(), 
                                    "3) The metadata/sample information (.txt or .csv) should contain variables for the subjects about host phenotypes, medical interventions, 
                                    disease status or environmental/behavioral factors, where rows are subjects and columns are variables (row names are subject IDs, and 
                                    column names are variable names).", br(), 
                                    "4) The phylogenetic tree (.tre or .nwk) should be a rooted tree. Otherwise, MiSurv automatically roots the tree through midpoint 
                                    rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.", br(), br(), 
                                    "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                                    The subjects should be matched and identical between feature table and metadata/sample information. 
                                    MiSurv will analyze only the matched features and subjects.", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT2 = p("You can download example microbiome data 'biom.zip'. This zip file contains four necessary data, feature table (otu.tab.txt), 
                                     taxonomic table (tax.tab.txt), metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre).", br(), 
                                     p(strong("Data description:"), "This example data are the public gut microbiome data (Zhang et al. 2018) we used in our real data applications (Gu et al., in review). The raw sequence data are deposited in QIITA (", tags$a(href = "https://qiita.ucsd.edu", "https://qiita.ucsd.edu"), ") with the ID number 10508 (", tags$a(href = "https://qiita.ucsd.edu/study/description/10508", "https://qiita.ucsd.edu/study/description/10508"), "). More detailed sample extraction and raw sequence data processing procedures are described in (Zhang et al. 2018; Gu et al., in review). "), 
                                     "> setwd('/yourdatadirectory/')", br(), br(), 
                                     "> otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE) ", br(), 
                                     "> tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)", br(), 
                                     " > sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE) ", br(),
                                     "> tree <- read.tree(file = 'tree.tre')", br(), br(), 
                                     " > sam.dat$T1Dweek",   HTML('&emsp;'),HTML('&nbsp;'),   "# Survival Time", br(), 
                                     " > sam.dat$T1D ", HTML('&emsp;'),HTML('&emsp;'),HTML('&emsp;'), "# Event", br(), 
                                     " > sam.dat$Antibiotics ", HTML('&emsp;'),  "# Treatment", br(), br(),
                                     "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, 
                                     and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                                     " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                                     " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                                     " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt")
  
  QC_KINGDOM_COMMENT = p("A microbial kingdom to be analyzed. Default is 'Bacteria' for 16S data. Alternatively, you can type 'Fungi' for ITS data 
                         or any other kingdom of interest for shotgun metagenomic data.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT1 = p("Remove subjects that have low library sizes (total read counts). Default is 3,000.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT2 = p("Library size: The total read count per subject.", style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT1 = p("Remove features (OTUs or ASVs) that have low mean relative abundances (Unit: %). Default is 0.002%.",style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT2 = p("Mean proportion: The average of relative abundances (i.e., proportions) per feature.", style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT1 = p('Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings. 
                            Multiple character strings should be separated by a comma. Default is "", "metagenome", "gut metagenome", "mouse gut metagenome".',
                            style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT2 = p('Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings (i.e., taxonomic names that contain 
                            the specified character strings). Multiple character strings should be separated by a comma. Default is "uncultured", "incertae", "Incertae",
                            "unidentified", "unclassified", "unknown".',style = "font-size:11pt")
  
  ALPHA_COMMENT = p("Calculate alpha-diversity indices: Richness (Observed), Shannon (Shannon, 1948), Simpson (Simpson, 1949), Inverse Simpson (Simpson, 1949), 
                    Fisher (Fisher et al., 1943), Chao1 (Chao, 1984), ACE (Chao and Lee, 1992), ICE (Lee and Chao, 1994), PD (Faith, 1992).")
  ALPHA_REFERENCES = p("1. Chao A, Lee S. Estimating the number of classes via sample coverage. J Am Stat Assoc. 1992:87:210-217.", br(),
                       "2. Chao A. Non-parametric estimation of the number of classes in a population. Scand J Stat. 1984:11:265-270.", br(),
                       "3. Faith DP. Conservation evaluation and phylogenetic diversity. Biol Conserv. 1992:61:1-10.", br(),
                       "4. Fisher RA, Corbet AS, Williams CB. The relation between the number of species and the number of individuals 
                       in a random sample of an animal population. J Anim Ecol. 1943:12:42-58.", br(),
                       "5. Lee S, Chao A. Estimating population size via sample coverage for closed capture-recapture models. Biometrics. 1994:50:1:88-97.", br(),
                       "6. Shannon CE. A mathematical theory of communication. Bell Syst Tech J. 1948:27:379-423 & 623-656.", br(),
                       "7. Simpson EH. Measurement of diversity. Nature 1949:163:688.", br())
  BETA_COMMENT = p("Calculate beta-diversity indices: Jaccard dissimilarity (Jaccard, 1912), Bray-Curtis dissimilarity (Bray and Curtis, 1957), Unweighted UniFrac distance 
                   (Lozupone and Knight, 2005), Generalized UniFrac distance (Chen et al., 2012), Weighted UniFrac distance (Lozupone et al., 2007).")
  BETA_REFERENCES = p("1. Bray JR, Curtis JT. An ordination of the upland forest communities of Southern Wisconsin. Ecol Monogr. 1957;27(32549).", br(),
                      "2. Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD., et al. Associating microbiome composition with environmental 
                      covariates using generalized UniFrac distances. Bioinformatics. 2012;28(16):2106-13.", br(),
                      "3. Jaccard P. The distribution of the flora in the alpine zone. New Phytol. 1912;11(2):37-50.", br(),
                      "4. Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and qualitative β-diversity measures lead to 
                      different insights into factors that structure microbial communities. Appl Environ Microbiol. 2007;73(5):1576-85.", br(),
                      "5. Lozupone CA, Knight R. UniFrac: A new phylogenetic method for comparing microbial communities. Appl Environ Microbiol. 2005;71(12):8228-35.")
  DATA_TRANSFORM_COMMENT = p("Transform the count (original) data into four different formats 1) count (rarefied) (Sanders, 1968), 2) proportion, 3) centered log ratio (CLR) (Aitchison, 1982), 4) arcsine-root for each taxonomic rank (phylum, class, order, family, genus, species).")
  DATA_TRANSFORM_REFERENCE = p("1. Aitchison J. The statistical analysis of compositional data. J R Statist Soc B. 1982;44:2:139-77")
  
  #Surv_Treat_Comment = p("Used to compare two categories of the selected binary variable.",style = "font-size:11pt")
  Surv_Time__Comment = p( "A variable for follow-up time (required).", style = "font-size:11pt")
  Surv_CensorComment = p("A binary indicator for censored (0) or event (1) (required).", style = "font-size:11pt")
  Surv_FollowComment = p("You can adjust the follow-up period (start time and end time) of interest (or leave it as it is to survey the entier follow-up period just as given in the data), where the censored/event indicator is also automatically adjusted to be suited to the selected follow-up period (optional).", style = "font-size:11pt")
  Follow_Time_Comment = p("You can reset the follow-up period of interest or leave it as it is.", style = "font-size:11pt")
  Num_PredComment = p("The number of taxa at each taxonomic rank to be displayed on the importance plot (default: 10).", style = "font-size:11pt")
}

# UI
{
  ui = dashboardPage(
    title = "MiSurv",
    dashboardHeader(title = span(TITLE, style = "float:left;font-size: 20px"), titleWidth = "100%"),
    dashboardSidebar(
      tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
      sidebarMenu(id = "side_menu",
                  menuItem("Home", tabName = "home", icon = icon("home")),
                  menuItem("Data Processing",  icon = icon("file-text-o"),
                           menuSubItem("Data Input", tabName = "step1", icon = icon("mouse")),
                           menuSubItem("Quality Control", tabName = "step2", icon = icon("chart-bar")),
                           menuSubItem("Data Transformation", tabName = "step3", icon = icon("calculator"))), #or, icon th-large
                  menuItem("Module 1", tabName = "SurvAnalysis", icon = icon("bar-chart-o")
                  ),
                  menuItem("Module 2",  icon = icon("chart-pie"),
                           menuSubItem("Alpha Diversity", tabName = "alphaDivanalysis", icon = icon("font")),
                           menuSubItem("Beta Diversity", tabName = "betaDivanalysis", icon = icon("bold")),
                           menuSubItem("Taxonomic Abundance", tabName = "taxaAnalysis", icon = icon("align-left"))),
                  menuItem("Module 3",  icon = icon("disease"),
                           menuSubItem("Alpha Diversity", tabName = "alphaDivanalysisSurv", icon = icon("font")),
                           menuSubItem("Beta Diversity", tabName = "betaDivanalysisSurv", icon = icon("bold")),
                           menuSubItem("Taxonomic Abundance", tabName = "taxaAnalysisSurv", icon = icon("align-left"))),
                  menuItem("Module 4", tabName = "RandomForest", icon = icon("tree")))),
    dashboardBody(
      tags$head(tags$style(HTML(".content { padding-top: 2px;}"))),
      tags$script(src = "fileInput_text.js"),
      useShinyjs(),
      shinyDashboardThemes(theme = "poor_mans_flatly"),
      uiOutput("themes"),
      tabItems(
        
        ##### HOME ####
        tabItem(tabName = "home",
                div(id = "homepage", br(), HOME_COMMENT, 
                    p(" ", style = "margin-bottom: 10px;"),
                    div(tags$img(src='MiSurv_workflow.png', height = 656, width = 1200), style="text-align: center;"),
                    br(),
                    tags$ul(
                      tags$li(HOME_COMMENT1), tags$li(HOME_COMMENT2), tags$li(HOME_COMMENT3), tags$li(HOME_COMMENT4), tags$li(HOME_COMMENT5),
                      style = "font-size:13pt"),
                    br(),
                    HOME_COMMENT6,
                    HOME_COMMENT7,
                    HOME_COMMENT8,
                    br()
                )),
        
        ##### DATA INPUT ####
        tabItem(tabName = "step1", br(),
                fluidRow(column(width = 6, style='padding-left:+20px',
                                box(
                                  width = NULL, status = "primary", solidHeader = TRUE,
                                  title = strong("Data Input"),
                                  selectInput("inputOption", h4(strong("Data Type?")), c("Choose one" = "", "Phyloseq", "Individual Data"), width = '30%'),
                                  div(id = "optionsInfo", tags$p("You can choose phyloseq or individual data.", style = "font-size:11pt"), style = "margin-top: -15px"),
                                  uiOutput("moreOptions"))),
                         column(width = 6, style='padding-left:0px', 
                                uiOutput("addDownloadinfo"), uiOutput("ref_micloud")))
        ),
        
        ##### QC ####
        tabItem(tabName = "step2", 
                br(), 
                fluidRow(column(width = 3, style = "padding-left:+20px", 
                                box(
                                  width = NULL, status = "primary", solidHeader = TRUE, 
                                  title = strong("Survival Data"), 
                                  h4(strong("Survival Time?", style = "color:black")),
                                  Surv_Time__Comment,
                                  p(" ", style = "margin-bottom: -20px;"),
                                  selectInput("surv.Time.select", label = "",
                                              c("Choose one" = "", ""), selected = FALSE, multiple = FALSE, selectize = FALSE, width = '70%'), 
                                  
                                  h4(strong("Censored/Event?", style = "color:black")),
                                  Surv_CensorComment,
                                  p(" ", style = "margin-bottom: -20px;"),
                                  selectInput("censor.select", label = "",
                                              c("Choose one" = "", ""), selected = FALSE, multiple = FALSE, selectize = FALSE, width = '70%'),
                                  
                                  p(" ", style = "margin-bottom: +20px;"),
                                  h4(strong("Adjust Follow-up Period?")),
                                  Surv_FollowComment,
                                  p(" ", style = "margin-bottom: -30px;"),
                                  sliderInput("follow.time", label = "", min=0, max=20, value = c(0,15), step = 1)
                                ),
                                box(
                                  width = NULL, status = "primary", solidHeader = TRUE, 
                                  title = strong("Microbiome Data"), 
                                  textInput("kingdom", h4(strong("Bacteria?")), value = "Bacteria"),
                                  QC_KINGDOM_COMMENT,
                                  tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                                  tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'),
                                  
                                  
                                  sliderInput("slider1", h4(strong("Library Size?")), min=0, max=10000, value = 3000, step = 1000),
                                  QC_LIBRARY_SIZE_COMMENT1,
                                  QC_LIBRARY_SIZE_COMMENT2,
                                  
                                  sliderInput("slider2", h4(strong("Mean Proportion?")), min = 0, max = 0.1, value = 0.002, step = 0.001,  post  = " %"),
                                  QC_MEAN_PROP_COMMENT1,
                                  QC_MEAN_PROP_COMMENT2,
                                  
                                  h4(strong("Erroneous Taxonomic Names?")),
                                  textInput("rem.str", label = "Complete Match", value = ""),
                                  QC_TAXA_NAME_COMMENT1,
                                  
                                  textInput("part.rem.str", label = "Partial Match", value = ""),
                                  QC_TAXA_NAME_COMMENT2,
                                  
                                  actionButton("run", (strong("Run!")), class = "btn-info"), br(), br(),
                                  uiOutput("moreControls")
                                )),
                         column(width = 9, style = "padding-left:+20px",
                                box(width = NULL, status = NULL, solidHeader = FALSE, 
                                    fluidRow(width = 12,
                                             status = "primary", solidHeader = TRUE, 
                                             valueBoxOutput("sample_Size", width = 3),
                                             valueBoxOutput("OTUs_Size", width = 3),
                                             valueBoxOutput("phyla", width = 3),
                                             valueBoxOutput("classes", width = 3)),
                                    fluidRow(width = 12, 
                                             status = "primary", solidHeader = TRUE,  #여깅
                                             valueBoxOutput("orders", width = 3),
                                             valueBoxOutput("families", width = 3),
                                             valueBoxOutput("genera", width = 3),
                                             valueBoxOutput("species", width = 3)),
                                    fluidRow(style = "position:relative",
                                             tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                                                    tabPanel("Histogram",
                                                             plotlyOutput("hist"),
                                                             sliderInput("binwidth", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                             chooseSliderSkin("Round", color = "#112446")),
                                                    tabPanel("Box Plot", 
                                                             plotlyOutput("boxplot"))),
                                             tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                                                    tabPanel("Histogram",
                                                             plotlyOutput("hist2"),
                                                             sliderInput("binwidth2", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                             chooseSliderSkin("Round", color = "#112446")),
                                                    tabPanel("Box Plot",
                                                             plotlyOutput("boxplot2"))))
                                ))),
                
                # mainPanel(width = 9,
                #           fluidRow(width = 12,
                #                    status = "primary", solidHeader = TRUE, 
                #                    valueBoxOutput("sample_Size", width = 3),
                #                    valueBoxOutput("OTUs_Size", width = 3),
                #                    valueBoxOutput("phyla", width = 3),
                #                    valueBoxOutput("classes", width = 3)),
                #           fluidRow(width = 12, 
                #                    status = "primary", solidHeader = TRUE,
                #                    valueBoxOutput("orders", width = 3),
                #                    valueBoxOutput("families", width = 3),
                #                    valueBoxOutput("genera", width = 3),
                #                    valueBoxOutput("species", width = 3)),
                #           fluidRow(style = "position:relative",
                #                    tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                #                           tabPanel("Histogram",
                #                                    plotlyOutput("hist"),
                #                                    sliderInput("binwidth", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                #                                    chooseSliderSkin("Round", color = "#112446")),
                #                           tabPanel("Box Plot", 
                #                                    plotlyOutput("boxplot"))),
                #                    tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                #                           tabPanel("Histogram",
                #                                    plotlyOutput("hist2"),
                #                                    sliderInput("binwidth2", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                #                                    chooseSliderSkin("Round", color = "#112446")),
                #                           tabPanel("Box Plot",
                #                                    plotlyOutput("boxplot2")))))
                
        ),
        
        ##### Data Transformation ####
        tabItem(tabName = "step3", br(),
                fluidRow(column(width = 6, style = 'padding-left:+20px',
                                box(title = strong("Data Transformation"), width = NULL, status = "primary", solidHeader = TRUE,
                                    ALPHA_COMMENT, 
                                    BETA_COMMENT, 
                                    DATA_TRANSFORM_COMMENT,
                                    actionButton("divCalcRun", (strong("Run!")), class = "btn-info"), 
                                    p(" ", style = "margin-bottom: +10px;"),
                                    p(strong("Attention:"), "Once you changed any setting in the preceding Data Input or Quality Control, you have to click this Run button again to use any of the following data analytic modules."),
                                    ),
                                uiOutput("divCalcDownload")),
                         column(width = 6, style='padding-left:0px',
                                box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                                    p("Alpha Diversity", style = "font-size:12pt"),
                                    ALPHA_REFERENCES,
                                    p("Beta Diversity", style = "font-size:12pt"),
                                    BETA_REFERENCES,
                                    p("Taxonomic Abundance", style = "font-size:12pt"),
                                    DATA_TRANSFORM_REFERENCE)))
                
        ),
        
        ##### Model 1: Survival Analysis ####
        tabItem(tabName = "SurvAnalysis", br(),
                sidebarLayout(
                  position = "left", 
                  div(style="width: 98%;",
                      sidebarPanel(width = 3,
                                   uiOutput("primvarsSurv"),
                                   uiOutput("prim_vars_types.Surv"),
                                   uiOutput("surviveTime"),
                                   uiOutput("censoring"),
                                   uiOutput("subgroup.sel"),
                                   uiOutput("covariates.Surv"), br(),
                                   uiOutput("referencesM1"))),
                  mainPanel(width = 9,
                            fluidRow(width = 8,
                                     status = "primary", solidHeader = TRUE, 
                                     div(style="margin-left: -1.2%; margin-right:-1.2%", valueBoxOutput("control_size", width = 3 ),
                                         valueBoxOutput("treatment_size", width = 3 ),
                                         valueBoxOutput("censored_size", width = 3 ),
                                         valueBoxOutput("uncensored_size", width = 3 ))),
                            #valueBoxOutput("significant_test", width = 3)),
                            fluidRow(
                              #div(style='width: 95%',
                              tabPanel(width = 8,"Survival Curve", #width = 12,
                                       uiOutput("surv_display_results")))))),#),
        
        
        ##### Model 2: ALPHA DIVERSITY ####
        tabItem(tabName = "alphaDivanalysis", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("primvars"),
                                                       uiOutput("prim_vars_types"),
                                                       uiOutput("covariates"), 
                                                       p(" ", style = "margin-bottom: +25px;"),
                                                       actionButton("runbtn_bin", (strong("Run!")), class = "btn-info"),
                                                       #br(), 
                                                       uiOutput("alpha_downloadTable"),
                                                       uiOutput("referencesM2.alpha"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("alpha_display_results"))))),
        
        ##### Model 2: BETA DIVERSITY ####
        tabItem(tabName = "betaDivanalysis", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("beta_primvar_cross"),
                                                       uiOutput("beta_prim_vars_types_cross"),
                                                       uiOutput("beta_covariates_cross"), 
                                                       actionButton("beta_runbtn_cross_bin", (strong("Run!")), class = "btn-info"),
                                                       #br(), 
                                                       uiOutput("beta_downloadTable"),
                                                       uiOutput("referencesM2.beta"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("beta_display_results_cross"))))),
        
        ##### Model 2: Taxa Analysis ####
        tabItem(tabName = "taxaAnalysis", br(),
                sidebarLayout( 
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("primvars_taxa"),
                                                       uiOutput("morePrimvar_opt_taxa"),
                                                       uiOutput("covariates_taxa"), 
                                                       actionButton("taxa_runbtn_bin", (strong("Run!")), class = "btn-info"),
                                                       #br(),
                                                       uiOutput("downloadTable_taxa"),
                                                       uiOutput("referencesM2.taxa"))),
                  
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     div(style='height:800px;overflow-y: scroll;', uiOutput("taxa_display_results")), br(),br(),
                                     
                                     uiOutput("taxa_display_dend"))
                            #column(9,  uiOutput("taxa_display_dend")),
                            #column(2, uiOutput("M2sig_taxon")))
                  ))),
        
        
        ##### Model 3: ALPHA DIVERSITY ####
        tabItem(tabName = "alphaDivanalysisSurv", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("subgroupSurv.A"),
                                                       uiOutput("subgroup.sel.A"),
                                                       uiOutput("surviveTimeCoxA"),
                                                       uiOutput("covariatesCoxA"), 
                                                       actionButton("runbtn_CoxA", (strong("Run!")), class = "btn-info"),
                                                       #br(),
                                                       uiOutput("alpha_surv_downloadTable"),
                                                       uiOutput("referencesM3.alpha"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("alpha_surv_display_resultsCoxA"))))),
        
        ##### Model 3: BETA  DIVERSITY ####
        tabItem(tabName = "betaDivanalysisSurv", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("subgroupSurv.B"),
                                                       uiOutput("subgroup.sel.B"),
                                                       uiOutput("surviveTimeCoxB"),
                                                       uiOutput("covariatesCoxB"), 
                                                       actionButton("runbtn_CoxB", (strong("Run!")), class = "btn-info"),
                                                       #br(),
                                                       uiOutput("beta_surv_downloadTable"),
                                                       uiOutput("referencesM3.beta"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     uiOutput("betaS_display_results_cross"))))),
        ##### Model 3: Taxonomic Abundance ####
        tabItem(tabName = "taxaAnalysisSurv", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("censor.selectCoxT"),
                                                       uiOutput("subgroupSurv.T"),
                                                       uiOutput("surviveTimeCoxT"),
                                                       uiOutput("covariatesCoxT"), 
                                                       actionButton("runbtn_CoxT", (strong("Run!")), class = "btn-info"),
                                                       #br(),
                                                       uiOutput("taxa_surv_downloadTable"),
                                                       uiOutput("referencesM3.taxa"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12, 
                                     div(style='height:800px;overflow-y: scroll;', uiOutput("surv_taxa_results")), br(),br(), #s_taxa_display_result 
                                     
                                     uiOutput("survival_taxa_display_dend"))
                            
                  ))),
        
        ##### Model 4: Prediction ####
        tabItem(tabName = "RandomForest", br(),
                sidebarLayout(
                  position = "left",
                  div(style="width: 98%;",sidebarPanel(width = 3,
                                                       uiOutput("model4_data_format"),
                                                       uiOutput("subgroupSurv.4"),
                                                       uiOutput("method.4"),
                                                       actionButton("runbtn_model4", (strong("Run!")), class = "btn-info"),
                                                       #br(),
                                                       uiOutput("downloadTable_m4"),
                                                       uiOutput("referencesM4"))),
                  mainPanel(width = 9,
                            fluidRow(width = 12,
                                     div(style='height:800px;overflow-y: scroll;', uiOutput("surv4_display_results"))))))
      )
    )
  )
}

# Server
server = function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  ## load example data ####
  env <- new.env()
  nm <- load(file = "Data/biom.Rdata", env)[1]
  biom <- env[[nm]]
  
  ori.biom <- biom
  otu.tab <- otu_table(ori.biom)
  tax.tab <- tax_table(ori.biom)
  tree <- phy_tree(ori.biom)
  sam.dat <- sample_data(ori.biom)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("biom",".Rdata", sep = "")
    },
    content = function(file1) {
      save(biom, file = file1)
    })
  
  output$downloadZip <- downloadHandler(
    filename = function() {
      paste("biom",".zip", sep = "")
    },
    content <- function(fname) {
      temp <- setwd(tempdir())
      on.exit(setwd(temp))
      dataFiles = c("otu.tab.txt", "tax.tab.txt", "sam.dat.txt" ,"tree.tre")
      write.table(otu.tab, "otu.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(tax.tab, "tax.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(sam.dat, "sam.dat.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.tree(tree, "tree.tre")
      zip(zipfile=fname, files=dataFiles)
    })
  
  ## variable define ####
  infile = reactiveValues(biom = NULL, qc_biom = NULL, qc_biomNA = NULL, rare_biom = NULL, rare_biomNA = NULL, is.mon = NULL, ori.var = NULL)
  ds.Ks <- reactiveValues(res = NULL)
  chooseData = reactiveValues(sam.dat = NULL, mon.sin.rev.bin.con = NULL, prim_vars = NULL, alpha.div = NULL,
                              alpha.div.rare = NULL, alpha.div.qc = NULL, taxa.out = NULL, taxa.outNA = NULL, tax.tab = NULL, tax.tabNA = NULL, 
                              prim_varsSurv = NULL, surv.time = NULL, nt.selected.cox = NULL,
                              prim_varsCox = NULL, nt.selected.coxA = NULL, nt.selected.sub = NULL, nt.selected.cov = NULL,
                              nt.selected.bin = NULL, nt.selected.con = NULL, #nt.selected.prim = NULL, trigger = NULL, 
                              sam_dat = NULL, surv.Time.select = NULL, censor.select = NULL, follow.time = NULL,
                              subgroup = NULL, subgroup.sel = NULL, pick_subgroup = NULL, NAadded = NULL)
  
  is.results = reactiveValues(result = NULL)
  is.results.long = reactiveValues(result = NULL)
  multi.test = reactiveValues(boolval = FALSE)
  multi.test.long = reactiveValues(boolval = FALSE)
  alpha.categos <- reactiveValues(cat1 = NULL, cat2 = NULL)
  alpha.categos.long <- reactiveValues(cat1 = NULL, cat2 = NULL)
  alpha.data.results = reactiveValues(table.out = NULL, data.q.out = NULL, table.p.out = NULL)
  alpha.results = reactiveValues(bin.var = NULL, alpha_div = NULL, alpha.bin.sum.out = NULL)
  alpha.results.cont = reactiveValues(alpha.con.out = NULL, alpha.table.out = NULL)
  alpha.reg.results = reactiveValues(bin.var = NULL, cov.var = NULL, alpha.div = NULL)
  alpha.resultslong = reactiveValues(bin.var = NULL, id.var = NULL, alpha_div = NULL, alpha.bin.sum.out = NULL)
  alpha.noncovs_res = reactiveValues(con.var = NULL, id.var = NULL, alpha.div = NULL)
  data.alphaBin_res = reactiveValues(table.output = NULL, data.output = NULL, alpha.bin.sum.out = NULL, table.p_outbin = NULL)
  data.results.cont_long = reactiveValues(table.out = NULL, data.q.out = NULL, table_p.out = NULL, alpha.table.out = NULL)
  
  beta.data.results = reactiveValues(data.q.out = NULL, betaS.plot.info = NULL)
  betaS.data.results = reactiveValues(data.q.out = NULL, betaS.plot.info = NULL)
  beta.results = reactiveValues(result = NULL)
  beta.resultscont = reactiveValues(beta.cont.out = NULL)
  beta.data.results_long = reactiveValues(beta.bin.out = NULL)
  beta.resultscon_long = reactiveValues(beta.con.out = NULL)
  beta.categos <- reactiveValues(cat1 = NULL, cat2 = NULL)
  beta.categos.long <- reactiveValues(cat1 = NULL, cat2 = NULL)
  beta.down.results <- reactiveValues(CS = NULL, LONG = NULL)
  
  alphaS.down.results <- reactiveValues(alpha.cox.out = NULL, alpha.cox.forest.plot = NULL)
  betaS.down.results <- reactiveValues(CS = NULL, LONG = NULL, mirkatS.plot = NULL)
  
  
  taxa.categos <- reactiveValues(cat1 = NULL, cat2 = NULL)
  taxa.data.results = reactiveValues(data.q.out = NULL)
  taxa.results = reactiveValues(bin.var = NULL, cov.var = NULL, id.var = NULL, taxa = NULL, taxa.bin.sum.out = NULL, con.var = NULL, taxa.con.sum.out = NULL, lib.size = NULL)
  taxa.types = reactiveValues(dataType = NULL, dataType_Surv = NULL, regression = NULL, dataType_model4 = NULL)
  taxa.outputs = reactiveValues(DAoutput = NULL, DAoutput_or = NULL, DAoutputlong = NULL)
  
  rcol = reactiveValues(selected = "lightblue") 
  
  ## options to change theme ####
  observeEvent(input$selectTheme, {
    output$themes <- renderUI({
      if (input$selectTheme == "Flat Red") {
        shinyDashboardThemes(
          theme = "flat_red"
        )
      } else if (input$selectTheme == "Gray Dark") {
        shinyDashboardThemes(
          theme = "grey_dark"
        )
      } else if (input$selectTheme =="Gray Light") {
        shinyDashboardThemes(
          theme = "grey_light"
        )
      } else if (input$selectTheme =="Onenote") {
        shinyDashboardThemes(
          theme = "onenote"
        )
      } else if (input$selectTheme == "Poor Mans Flatly (Default)") {
        shinyDashboardThemes(
          theme = "poor_mans_flatly"
        )
      } else if (input$selectTheme == "Purple Gradient") {
        shinyDashboardThemes(
          theme = "purple_gradient"
        )
      }
    })
  }, ignoreInit = TRUE)
  ## options to input ####
  observeEvent(input$inputOption,{
    observe({
      if (input$inputOption == "Phyloseq") {
        shinyjs::hide(id = "optionsInfo")
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("phyloseqData", strong("Please upload your 'phyloseq' data (.Rdata, .rds)", style = "color:black"), 
                      accept = c(".Rdata", ".rds"), width = '80%'), div(style = "margin-top: -18px"),
            
            actionButton('Load_Phyloseq_Data', 'Upload', class = "btn-info"),
            br(),
            p(" ", style = "margin-top: 5px;"),
            p(strong("Attention: ", style = "color:black"), "You have to click this Upload button to perform following data processing and downstream data analyses."), 
            
            br(), 
            shinyjs::hidden(
              shiny::div(id = "phyloseqUpload_error",
                         shiny::tags$p("Please upload a Rdata file!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_PHYLOSEQ_COMMENT1
          )
        })
        
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data"), width = NULL, status = "primary", solidHeader = TRUE,
                downloadButton("downloadData", "Example Data", width = '30%', style = "background-color: red2"),
                br(),br(),
                INPUT_PHYLOSEQ_COMMENT2
            )
          )
        })
        
        output$ref_micloud <- renderUI({
          tagList(
            box(title = strong("Reference"), width = NULL, status = "primary", solidHeader = TRUE, 
                p("Zhang XS, Li J, Krautkramer KA, Badri M, Battaglia T, Borbet TC et al. Antibiotic-induced acceleration of type 1 diabetes alters maturation of innate intestinal immunity. eLife 2018:7:e37816."), 
                p("Gu W, Koh H, Jang HJ, Lee B, Kang B. MiSurv: An integrative web cloud platform for user-friendly microbiome data analysis with survival responses. (in review)")
            )
          )
        })
        
      } else if (input$inputOption == "Individual Data") {
        shinyjs::hide(id = "optionsInfo")
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("otuTable", strong("Please upload your feature (OTU or ASV) table (.txt, .csv, .biom)", style = "color:black"), 
                      accept = c(".txt", ".csv", ".biom"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("taxTable", strong("Please upload your taxonomic table (.txt, .tsv)", style = "color:black"), 
                      accept = c(".txt", ".tsv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("samData", strong("Please upload your metadata/sample information (.txt, .csv)", style = "color:black"), 
                      accept = c(".txt", ".csv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("tree", strong("Please upload your phylogenetic tree (.tre, .nwk)", style = "color:black"), 
                      accept = c(".tre", ".nwk"), width = '80%'), div(style = "margin-top: -15px"),
            
            
            actionButton('Load_Individual_Data', 'Upload', class = "btn-info"), 
            br(),
            p(" ", style = "margin-top: 5px;"),
            p(strong("Attention: ", style = "color:black"), "You have to click this Upload button to perform following data processing and downstream data analyses."), 
            
            br(), 
            shinyjs::hidden(
              shiny::div(id = "textfilesUpload_error",
                         shiny::tags$p("Please upload txt and tre files!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_INDIVIDUAL_DATA_COMMENT
          )
        })
        
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data"), width = NULL, status = "primary", solidHeader = TRUE,
                downloadButton("downloadZip", "Example Data", width = '30%', style = "background-color: red2"),
                br(),br(),
                INPUT_INDIVIDUAL_DATA_COMMENT2
            )
          )
        })
        
        
        output$ref_micloud <- renderUI({
          tagList(
            box(title = strong("Reference"), width = NULL, status = "primary", solidHeader = TRUE, 
                p("Zhang XS, Li J, Krautkramer KA, Badri M, Battaglia T, Borbet TC et al. Antibiotic-induced acceleration of type 1 diabetes alters maturation of innate intestinal immunity. eLife 2018:7:e37816."), 
                p("Gu W, Koh H, Jang HJ, Lee B, Kang B. MiSurv: An integrative web cloud platform for user-friendly microbiome data analysis with survival responses. (in review)")
            )
          )
        })
      }
    })
    
  }, ignoreInit = TRUE, once = TRUE, ignoreNULL = TRUE)
  
  observe({
    toggleState("Load_Phyloseq_Data", !is.null(input$phyloseqData))
    toggleState("Load_Individual_Data", 
                !(is.null(input$otuTable) | is.null(input$taxTable) | is.null(input$samData) | is.null(input$tree)))
    toggleState("run", !(is.null(infile$biom)))
    toggleState("skip", !is.null(infile$biom))
    toggleState("slider1", !is.null(infile$biom))
    toggleState("slider2", !is.null(infile$biom))
    toggleState("kingdom", !is.null(infile$biom))
    toggleState("binwidth", !is.null(infile$biom))
    toggleState("binwidth2", !is.null(infile$biom))
    toggleState("surv.Time.select", !is.null(infile$biom))
    toggleState("follow.time", !is.null(infile$biom))
    toggleState("subgroup", !is.null(infile$biom))
    toggleState("subgroup.sel", !is.null(infile$biom))
    toggleState("censor.select", !is.null(infile$biom))
    
    toggleState("divCalcRun", !is.null(infile$rare_biom))
    #toggleState("datTransRun", !is.null(infile$rare_biom))
    
    toggleState("surv.KMrun", !is.null(infile$rare_biom))
  })
  
  observeEvent(input$Load_Phyloseq_Data, {
    
    if (!is.null(input$phyloseqData)) {
      dataInfile  = reactive({
        phyloseq.data = input$phyloseqData
        ext <- tools::file_ext(phyloseq.data$datapath)
        
        req(phyloseq.data)
        if (ext == "Rdata") {
          phyloseq.dataPath = phyloseq.data$datapath
          e = new.env()
          name <- load(phyloseq.dataPath, envir = e)
          data <- e[[name]]
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          
          colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } else if (ext == "rds") {
          phyloseq.dataPath = phyloseq.data$datapath
          data <- readRDS(phyloseq.dataPath)
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } else {
          shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade")
          shinyjs::delay(5000, shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade"))
          return(NULL)
        }
      })
    } else {
      return(NULL)
    }
    
    if (is.null(dataInfile)) {
      infile$biom <- NULL
      infile$qc_biom <- NULL
      infile$rare_biom = NULL
    } else {
      infile$biom <- dataInfile()
      infile$qc_biom <- dataInfile()
      infile$rare_biom = NULL
    }
    
    updateTabsetPanel(session, "side_menu",
                      selected = "step2")
    rcol$selected = "lightblue"
    
    if (!is.null(infile$biom)) QC$resume()
  })
  observeEvent(input$Load_Individual_Data, {
    shinyjs::disable("Load_Individual_Data")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(3/10, message = "File Check")
        if (!is.null(input$otuTable) & !is.null(input$taxTable) & !is.null(input$samData) & !is.null(input$tree)) {
          dataInfile  = reactive({
            otu.table = input$otuTable
            ext1 <- tools::file_ext(otu.table$datapath)
            
            tax.table = input$taxTable
            ext2 <- tools::file_ext(tax.table$datapath)
            
            sam.data = input$samData
            ext3 <- tools::file_ext(sam.data$datapath)
            
            tree.data = input$tree
            ext4 <- tools::file_ext(tree.data$datapath)
            
            req(otu.table, tax.table, sam.data, tree.data)
            if ((ext1 == "txt"| ext1 == "csv" | ext1 == "biom") & (ext2 == "txt" | ext2 == "tsv") &
                (ext3 == "txt" | ext3 == "csv") & (ext4 == "tre" | ext4 == "nwk")) {
              otu.table.path = otu.table$datapath
              tax.table.path = tax.table$datapath
              sam.data.path = sam.data$datapath
              tree.data.path = tree.data$datapath
              
              if (ext1 == "txt") {
                otu.tab <- read.table(otu.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext1 == "csv") {
                otu.tab <- read.csv(otu.table.path, check.names = FALSE)
                rownames(otu.tab) = otu.tab[,1];otu.tab = otu.tab[,-1]
              } else if (ext1 == "biom") {
                biom <- read_biom(otu.table.path)
                otu.tab <- as.matrix(biom_data(biom))
              }
              
              if (ext2 == "txt") {
                tax.tab <- read.table(tax.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext2 == "tsv") {
                tax.tab <- read.table(tax.table.path, header=TRUE, sep="\t")
                tax.tab = preprocess.tax.tab(tax.tab)
              }
              
              if (ext3 == "txt") {
                sam.dat <- read.table(sam.data.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext3 == "csv") {
                sam.dat <- read.csv(sam.data.path, check.names = FALSE)
                rownames(sam.dat) = sam.dat[,1]
                sam.dat = sam.dat[,-1]
              }
              
              if (ext4 == "tre") {
                tree <- read.tree(file = tree.data.path)
              } else if (ext4 == "nwk") {
                tree <- read.tree(file = tree.data.path) 
              }
              
              otu.tab <- otu_table(otu.tab, taxa_are_rows = TRUE)
              tax.tab <- tax_table(as.matrix(tax.tab))
              sam.dat <- sample_data(sam.dat)
              tree <- phy_tree(tree)
              
              if (sum(colnames(otu.tab) %in% rownames(sam.dat)) < sum(rownames(otu.tab) %in% rownames(sam.dat))) {
                otu.tab = t(otu.tab)
              }
              
              incProgress(3/10, message = "Validating")
              validate(
                if (biom.check.samples(otu.tab, sam.dat)) {
                  if (biom.check.otu(otu.tab, tax.tab, tree)) {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data. And
                                        there is no common OTUs among OTU/feature table, taxonomic table and tree tip labels"),
                                     type = "error")
                  } else {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data"),
                                     type = "error")
                  }
                } else if (biom.check.otu(otu.tab, tax.tab, tree)) {
                  showNotification(h4("Error: There is no common OTUs among OTU/feature table, taxonomic table and tree tip labels"),
                                   type = "error")
                } else {
                  NULL
                }
              )
              
              incProgress(1/10, message = "Merging")
              biomData <- merge_phyloseq(otu.tab, tax.tab, sam.dat, tree)
              return(biomData)
            } else {
              shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(5000, shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade"))
              return(NULL)
            }
          })
        } else {
          return(NULL)
        }
        
        if (is.null(dataInfile)) {
          infile$biom <- NULL
          infile$qc_biom <- NULL
          infile$rare_biom = NULL
        } else {
          infile$biom <- dataInfile()
          infile$qc_biom <- dataInfile()
          infile$rare_biom = NULL
        }
        
        updateTabsetPanel(session, "side_menu",
                          selected = "step2")
        rcol$selected = "lightblue"
        
        if (!is.null(infile$biom)) QC$resume()
      })
    shinyjs::enable("Load_Individual_Data")
  })
  
  ######################################
  # Quality control and transformation #
  ######################################
  # This reactive expression stores the input data from either the individual data or phyloseq data
  QC = observe(suspended = T,{
    taxa.results$lib.size <- lib.size.func(infile$biom)$lib.size
    
    # Plots graphs using example infile data
    output$hist <- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      plot_ly(x = ~lib_size, nbinsx = input$binwidth,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Library Size", zeroline = FALSE))
    })
    
    output$hist2 <- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      plot_ly(x = ~mean_prop, nbinsx = input$binwidth2,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Mean Proportion", zeroline = FALSE))
    })
    
    output$boxplot<- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      
      plot_ly(x = ~lib_size, type = "box", notched=TRUE, name = "Library Size",
              color = ~"lib_size", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    output$boxplot2<- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      
      plot_ly(x = ~mean_prop, type = "box", notched=TRUE, name = "Mean Proportion",
              color = ~"mean_prop", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    ## Number of Taxonomic Rank for biom before QC
    num_tax.rank = reactive({
      tax.tab = tax_table(infile$qc_biom)
      num.tax.rank(tax.tab)
    })
    
    ## Fills value boxes using example biom data
    output$sample_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.sams), style = "font-size: 75%;"),
        "Sample Size", icon = icon("user-circle"), color = "fuchsia")
    })
    
    output$OTUs_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.otus), style = "font-size: 75%;"),
        "Number of Features", icon = icon("dna"), color = "aqua")
    })
    
    output$phyla <- renderValueBox({
      num.phyla = num_tax.rank()[1]
      valueBox(
        value = tags$p(paste0(num.phyla), style = "font-size: 75%;"),
        "Number of Phyla", icon = icon("sitemap"), color = "orange")
    })
    
    output$classes <- renderValueBox({
      num.classes = num_tax.rank()[2]
      valueBox(
        value = tags$p(paste0(num.classes), style = "font-size: 75%;"),
        "Number of Classes", icon = icon("sitemap"), color = "purple")
    })
    
    output$orders <- renderValueBox({
      num.orders = num_tax.rank()[3]
      valueBox(
        value = tags$p(paste0(num.orders), style = "font-size: 75%;"),
        "Number of Orders", icon = icon("sitemap"), color = "blue")
    })
    
    output$families <- renderValueBox({
      num.families = num_tax.rank()[4]
      valueBox(
        value = tags$p(paste0(num.families), style = "font-size: 75%;"),
        "Number of Families", icon = icon("sitemap"), color = "red")
    })
    
    output$genera <- renderValueBox({
      num.genera = num_tax.rank()[5]
      valueBox(
        value = tags$p(paste0(num.genera), style = "font-size: 75%;"),
        "Number of Genera", icon = icon("sitemap"), color = "lime")
    })
    
    output$species <- renderValueBox({
      num.species = num_tax.rank()[6]
      valueBox(
        value = tags$p(paste0(num.species), style = "font-size: 75%;"),
        "Number of Species", icon = icon("sitemap"), color = "teal" )
    })
    
    ## This event handler checks whether there is an input file and updates the slider options accordingly
    maxi.slider1 = as.numeric(lib.size.func(infile$qc_biom)$lib.size.sum["3rd quartile"])
    max.mean.prop = as.numeric(mean.prop.func(infile$qc_biom)$mean.prop.sum["3rd quartile"])
    maxi.slider2 = round(max.mean.prop, digits = 6)
    
    if (maxi.slider2 < 2e-05) {
      maxi.slider2 = 2e-05
    }
    
    updateSliderInput(session, "slider1", min = 0, max = round(maxi.slider1,-3))
    updateSliderInput(session, "slider2", min = 0, max = maxi.slider2*100)
    
    observeEvent(c(input$surv.Time.select, input$censor.select),{
      if (input$surv.Time.select == c("Choose one" = "", "") || input$censor.select == c("Choose one" = "", "") ){
        disable("run")
      } else if (input$surv.Time.select != c("Choose one" = "", "") && input$censor.select != c("Choose one" = "", "") ){
        enable("run")
      }
    })
    
    observeEvent(input$surv.Time.select, {
      if(is.null(chooseData$follow.time)) {
        updateSliderInput(session, "follow.time", min = 0, max = max(sample_data(infile$biom)[[input$surv.Time.select]]), 
                          value = c(0, quantile(sample_data(infile$biom)[[input$surv.Time.select]])[[4]]))
      } else {
        updateSliderInput(session, "follow.time", min = 0, max = max(sample_data(infile$biom)[[input$surv.Time.select]]), 
                          value = as.numeric(chooseData$follow.time))
      }
    })
  })
  
  observeEvent(infile$biom, {
    infile$is.mon <- is.mon.sin.rev.bin.con(sample_data(infile$biom))
    infile$ori.var <- surv.pri.func(sample_data(infile$biom), infile$is.mon)
    # if(is.null(input$surv.Time.select) || nchar(input$surv.Time.select) == 0) {
    updateSelectInput(session, "surv.Time.select", choices = c("Choose one" = "", sort(infile$ori.var[[2]])), selected = "")
    updateSelectInput(session, "censor.select", label = "", c("Choose one" = "", sort(infile$ori.var[[1]])), selected = "")
    # }
  })
  
  ######################################
  ######################################
  #########  Data Analysis   ###########
  ######################################
  ######################################
  
  observeEvent(input$run, {
    observeEvent(chooseData$prim_varsSurv, {
      ######################################
      ######################################
      ######## Survival Analysis ###########
      ######################################
      ######################################
      
      if (length(chooseData$prim_varsSurv[[1]]) > 0){
        
        # 1. Treatment
        output$primvarsSurv <- renderUI({
          tagList(
            selectInput("primvar.Surv.select", label = h4(strong("Treatment Variable?", style = "color:black")),
                        c("Choose one" = "", sort(chooseData$nt.selected.bin)), selected = sort(chooseData$nt.selected.bin)[1], width = '80%')) 
        })
        
        output$prim_vars_types.Surv <- renderUI({
          tagList(
            uiOutput("morePrim_Surv"))
        })
        
        
      } else {
        showNotification("Error: No analysis available - Binary variable not found.", duration = 15,
                         type = "error")
      }
      
      observeEvent(input$primvar.Surv.select, {
        
        Surv.bin.cat.ref.ori.out <- alpha.bin.cat.ref.ori.func(chooseData$sam.dat, input$primvar.Surv.select)
        
        # 2. Treatment Name Change
        output$morePrim_Surv <- renderUI({
          tagList(
            h4(strong("Rename Categories?", style = "color:black")),
            p("You can rename the categories of treatment variable. MiSurv keeps up to 8 characters on graphs.", style = "font-size:11pt"),
            textInput("prim.var.rename1.Surv", label = (paste0("Reference: ",Surv.bin.cat.ref.ori.out[1])), value = Surv.bin.cat.ref.ori.out[1], width = '80%'),
            textInput("prim.var.rename2.Surv", label = (paste0("Comparison: ",Surv.bin.cat.ref.ori.out[2])), value = Surv.bin.cat.ref.ori.out[2], width = '80%'))
        })
        
        
        nt.selected.prim <- c(chooseData$nt.selected.bin, chooseData$nt.selected.con)[!c(chooseData$nt.selected.bin, chooseData$nt.selected.con) %in% input$primvar.Surv.select]
        
        output$covariates.Surv <- renderUI({
          tagList(
            p(" ", style = "margin-bottom: -5px;"),
            
            prettyRadioButtons("covariates.surv",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                               animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
            
            shinyjs::hidden(
              shiny::div(id = "covariates_variables.surv", style = "margin-left: 2%",
                         prettyCheckboxGroup("covariatesOptions.surv"," Please select covariate(s):", status = "primary",
                                             nt.selected.prim, width = '70%'))),
            
            uiOutput("chooseTest.surv"),
            
            uiOutput("chooseTest.surv2"),
            
            actionButton("surv.KMrun", (strong("Run!")), class = "btn-info")
          )
        })
        
        output$chooseTest.surv2 <- renderUI({
          tagList(
            p(" ", style = "margin-top: -20px;"),
            selectInput("sig.test.Surv", label = h4(strong("", style = "color:black")),
                        choices = c("Log-rank test", "Wilcoxon test"),
                        selected = "Log-rank test", width = '80%')
          )
        }) 
        
        if(!is.null(input$pick_KM_Cox)) {
          if(input$pick_KM_Cox == "Kaplan-Meier analysis (default)") {
            shinyjs::show("chooseTest.surv2")
          }
          else if(input$pick_KM_Cox == "Cox model"){
            shinyjs::hide("chooseTest.surv2")
          }
        }
      })
      
      observeEvent(input$covariates.surv,{
        if (input$covariates.surv == "None") {
          shinyjs::show("chooseTest.surv2")
          shinyjs::hide("covariates_variables.surv")
          
          output$chooseTest.surv <- renderUI({
            tagList(
              
              p(" ", style = "margin-top: 5px;"),
              
              h4(strong("Method?", style = "color:black")),
              
              p(" ", style = "margin-bottom: -20px;"),
              
              radioButtons("pick_KM_Cox", label = "", #label = h4(strong("Method?", style = "color:black")),
                           choices = c("Kaplan-Meier analysis (default)", "Cox model"),
                           selected = "Kaplan-Meier analysis (default)", width = '80%')
            )
          })
        }
        else if (input$covariates.surv == "Covariate(s)") {
          
          shinyjs::hide("chooseTest.surv2")
          
          shinyjs::show("covariates_variables.surv")
          
          output$chooseTest.surv <- renderUI({
            tagList(
              p(" ", style = "margin-bottom: 5px;"),
              
              h4(strong("Method?", style = "color:black")),
              
              p(" ", style = "margin-bottom: -20px;"),
              
              radioButtons("pick_KM_Cox", label = "", #h4(strong("Method?", style = "color:black")),
                           choices = c("Cox model"),
                           selected = "Cox model", width = '80%')
            )
          })
        }
      })
    })
  })
  
  observeEvent(input$pick_KM_Cox, {
    if(input$pick_KM_Cox == "Kaplan-Meier analysis (default)") {shinyjs::show("chooseTest.surv2")}
    else if(input$pick_KM_Cox == "Cox model") {shinyjs::hide("chooseTest.surv2")}
  })
  
  observeEvent(chooseData$alpha.div, {
    
    ######################################
    ######################################
    ######### Alpha diversity ############
    ######################################
    ######################################
    
    ######################################
    ######## Alpha Surv Analysis #########
    ######################################
    
    output$subgroupSurv.A <- renderUI({
      tagList(
        
        h4(strong("Subgroup Analysis?", style = "color:black")), 
        
        p("You can select a subgroup of interest using a category of a nominal variable (usually the treatment variable) to perform subgroup analysis. If you select a subgroup, only the subjects in that subgroup are retained in the following analysis. Otherwise, all the subjects given in the data are retained (optional).", style = "font-size:11pt"),
        
        prettyRadioButtons("subgroup.a", label = NULL, status = "primary", icon = icon("check"),
                           animation = "jelly", choices = c("No", "Yes"), selected="No", width = '70%'),
        p(" ", style = "margin-bottom: -25px;"),
        shinyjs::hidden(
          tagList(
            shiny::div(id = "subgroup.sel.a", style = "margin-left: 5%",
                       selectInput("subgroup.sel.a", label = "", 
                                   c("Choose one" = "", "null"), selected = "null", width = '70%'))
          )
        ),
        p(" ", style = "margin-bottom: -10px;"),
        uiOutput("pick_subgroup.a"),
        p(" ", style = "margin-bottom: +20px;"),
      )
    })
    
    observeEvent(input$subgroup.a,{
      if (input$subgroup.a == "Yes"){
        cat.for.subgroup =  chooseData$prim_varsCox[[1]][ ! chooseData$prim_varsCox[[1]]%in% input$censor.select ]
        
        if(!is.null(input$subgroup.sel.a)){
          if(input$subgroup.sel.a %in% cat.for.subgroup) {
            updateSelectInput(session, "subgroup.sel.a", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                              selected = input$subgroup.sel.a)
          } else {
            updateSelectInput(session, "subgroup.sel.a", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                              selected = sort(cat.for.subgroup)[1])
          } 
        } else {
          updateSelectInput(session, "subgroup.sel.a", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                            selected = sort(cat.for.subgroup)[1])
        }
        shinyjs::show("subgroup.sel.a")
        shinyjs::show("pick_subgroup.a")
      }
      else {
        shinyjs::hide("subgroup.sel.a")
        shinyjs::hide("pick_subgroup.a")
      }
    })
    
    
    observeEvent(input$subgroup.sel.a, {
      if(input$subgroup.sel.a %in% infile$ori.var[[1]]){
        
        
        subgroup.val <- alpha.bin.cat.ref.ori.func(chooseData$sam.dat, input$subgroup.sel.a)
        
        output$pick_subgroup.a <- renderUI({
          tagList(
            p(" ", style = "margin-bottom: -5px;"),
            shiny::div(id = "pick_subgroup.a", style = "margin-left: 5%",
                       radioButtons("pick_subgroup.a", label =  "Please select a subgroup:", 
                                    choices = subgroup.val,
                                    selected = subgroup.val[1], width = '70%')),
            p(" ", style = "margin-bottom: -20px;")
          )
        })
      }
    })
    
    output$covariatesCoxA <- renderUI({
      tagList(
        p(" ", style = "margin-top: +25px;"),
        prettyRadioButtons("covariatesCoxA",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                           animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
        
        shinyjs::hidden(
          shiny::div(id = "covariates_variablesCoxA", style = "margin-left: 2%",
                     prettyCheckboxGroup("covariatesOptionsCoxA","", status = "primary", "null", width = '70%'))),
        
        
        p(" ", style = "margin-top: 5px;"),
        
        h4(strong("Method?", style = "color:black")),
        
        p(" ", style = "margin-bottom: -20px;"),
        
        radioButtons("surv3.alpha.method", label = "", #h4(strong("Method?", style = "color:black")),
                     c("Cox model"), selected = "Cox model", width = '80%'),
        
        p(" ", style = "margin-bottom: -10px;")
        
      )
    })
    
    observeEvent(input$subgroup.sel.a,{
      
      options_bin <- chooseData$prim_varsSurv[[1]][ !(chooseData$prim_varsSurv[[1]] %in% c(input$censor.select, input$subgroup.sel.a))]
      options_con <- chooseData$prim_varsSurv[[2]][ !(chooseData$prim_varsSurv[[2]] %in% c(input$surv.Time.select, input$subgroup.sel.a))]
      
      updatePrettyCheckboxGroup(session, "covariatesOptionsCoxA", " Please select covariate(s)", sort(c(options_bin, options_con)), NULL)
    })
    
    
    observeEvent(input$covariatesCoxA,{
      if (input$covariatesCoxA == "Covariate(s)") {
        
        shinyjs::show("covariates_variablesCoxA")
        
      } else if (input$covariatesCoxA == "None") {
        
        shinyjs::hide("covariates_variablesCoxA")
        
      }
    })
    
    ######################################
    ######################################
    ########## Beta diversity ############
    ######################################
    ######################################
    
    ######################################
    ########## Beta Surv Analysis ########
    ######################################
    output$subgroupSurv.B <- renderUI({
      tagList(
        
        h4(strong("Subgroup Analysis?", style = "color:black")), 
        p("You can select a subgroup of interest using a category of a nominal variable (usually the treatment variable) to perform subgroup analysis. If you select a subgroup, only the subjects in that subgroup are retained in the following analysis. Otherwise, all the subjects given in the data are retained (optional).", style = "font-size:11pt"),
        
        prettyRadioButtons("subgroup.b", label = NULL, status = "primary", icon = icon("check"),
                           animation = "jelly", choices = c("No", "Yes"), selected="No", width = '70%'),
        p(" ", style = "margin-bottom: -25px;"),
        shinyjs::hidden(
          tagList(
            shiny::div(id = "subgroup.sel.b", style = "margin-left: 5%",
                       selectInput("subgroup.sel.b", label = "", 
                                   c("Choose one" = "", "null"), selected = "null", width = '70%'))
          )
        ),
        p(" ", style = "margin-bottom: -10px;"),
        uiOutput("pick_subgroup.b"),
        p(" ", style = "margin-bottom: +20px;"),
      )
    })
    
    observeEvent(input$subgroup.b,{
      if (input$subgroup.b == "Yes"){
        cat.for.subgroup = chooseData$prim_varsCox[[1]][ !chooseData$prim_varsCox[[1]]%in% input$censor.select ]
        
        if(!is.null(input$subgroup.sel.b)){
          if(input$subgroup.sel.b %in% cat.for.subgroup) {
            updateSelectInput(session, "subgroup.sel.b", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                              selected = input$subgroup.sel.b)
          } else {
            updateSelectInput(session, "subgroup.sel.b", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                              selected = sort(cat.for.subgroup)[1])
          } 
        } else {
          updateSelectInput(session, "subgroup.sel.b", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                            selected = sort(cat.for.subgroup)[1])
        }
        shinyjs::show("subgroup.sel.b")
        shinyjs::show("pick_subgroup.b")
      }
      else {
        shinyjs::hide("subgroup.sel.b")
        shinyjs::hide("pick_subgroup.b")
      }
    })
    
    
    observeEvent(input$subgroup.sel.b, {
      if(input$subgroup.sel.b %in% infile$ori.var[[1]]){
        
        subgroup.val <- alpha.bin.cat.ref.ori.func(chooseData$sam.dat, input$subgroup.sel.b)
        
        output$pick_subgroup.b <- renderUI({
          tagList(
            p(" ", style = "margin-bottom: -5px;"),
            shiny::div(id = "pick_subgroup.b", style = "margin-left: 5%",
                       radioButtons("pick_subgroup.b", label =  "Please select a subgroup:", 
                                    choices = subgroup.val,
                                    selected = subgroup.val[1], width = '70%')),
            p(" ", style = "margin-bottom: -20px;")
          )
        })
      }
    })
    
    output$covariatesCoxB<- renderUI({
      tagList(
        p(" ", style = "margin-top: +25px;"),
        prettyRadioButtons("covariatesCoxB",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                           animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
        
        shinyjs::hidden(
          shiny::div(id = "covariates_variablesCoxB", style = "margin-left: 2%",
                     prettyCheckboxGroup("covariatesOptionsCoxB","", status = "primary", "null", width = '70%'))),
        
        
        p(" ", style = "margin-top: 5px;"),
        
        h4(strong("Method?", style = "color:black")),
        
        p(" ", style = "margin-bottom: -20px;"),
        
        radioButtons("surv3.beta.method", label = "", c("MiRKAT-S"), selected = "MiRKAT-S", width = '80%'),
        
        p(" ", style = "margin-bottom: -10px;")
        
      )
    })
    
    observeEvent(input$subgroup.sel.b,{
      
      options_bin <- chooseData$prim_varsSurv[[1]][ !(chooseData$prim_varsSurv[[1]] %in% c(input$censor.select, input$subgroup.sel.b))]
      options_con <- chooseData$prim_varsSurv[[2]][ !(chooseData$prim_varsSurv[[2]] %in% c(input$surv.Time.select, input$subgroup.sel.b))]
      
      updatePrettyCheckboxGroup(session, "covariatesOptionsCoxB", " Please select covariate(s)", sort(c(options_bin, options_con)), NULL)
    })
    
    observeEvent(input$covariatesCoxB,{
      if (input$covariatesCoxB == "Covariate(s)") {
        
        shinyjs::show("covariates_variablesCoxB")
        
      } else if (input$covariatesCoxB == "None") {
        
        shinyjs::hide("covariates_variablesCoxB")
        
      }
    })
  })
  
  observeEvent(chooseData$alpha.div,{
    
    ######################################
    ######################################
    ######### Alpha diversity ############
    ######################################
    ######################################
    
    ######################################
    ### Cross-sectional Data Analysis ####
    ######################################
    
    output$primvars <- renderUI({
      tagList(
        selectInput("primvar", label = h4(strong("Treatment Variable?", style = "color:black")),
                    c("Choose one" = "", sort(chooseData$nt.selected.bin)), selected = sort(chooseData$nt.selected.bin)[1], width = '80%'))
    })
    
    
    output$prim_vars_types <- renderUI({
      tagList(
        uiOutput("morePrimvar_opt"))
    })
    
    observeEvent(input$primvar,{
      ## user selects whether to use rarefied or non rarefied biom data
      
      if (input$primvar %in% chooseData$nt.selected.bin) {
        is.results$result = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar)
        
        nt.selected_cat =  chooseData$prim_varsSurv[[1]][ chooseData$prim_varsSurv[[1]] != input$primvar]
        chooseData$nt.selected.sub = nt.selected_cat
        
        # if variable is binary, perform data analysis
        if (is.results$result == "Binary") {
          alpha.categos$cat1 = alpha.bin.cat.func(chooseData$sam.dat, input$primvar)[1]
          alpha.categos$cat2 = alpha.bin.cat.func(chooseData$sam.dat, input$primvar)[2]
          
          output$morePrimvar_opt <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")),
              p("You can rename the categories of primary variable. MiSurv keeps up to 8 characters on graphs.", style = "font-size:11pt"), 
              textInput("alphaCat1", label = (paste0("Reference: ",alpha.categos$cat1)), value = alpha.categos$cat1, width = '80%'),
              textInput("alphaCat2", label = (paste0("Comparison: ",alpha.categos$cat2)), value = alpha.categos$cat2, width = '80%'))
          })
          
          nt.selected.prim <- c(chooseData$nt.selected.bin, chooseData$nt.selected.con)[!c(chooseData$nt.selected.bin, chooseData$nt.selected.con) %in% input$primvar.Surv.select]
          
          output$covariates <- renderUI({
            tagList(
              p(" ", style = "margin-bottom: -5px;"),
              prettyRadioButtons("covariates",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                                 animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
              shinyjs::hidden(
                shiny::div(id = "covariates_variables", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariatesOptions"," Please select covariate(s)", status = "primary",
                                               nt.selected.prim, width = '70%'))),
              
              uiOutput("chooseTest"),
              
              p(" ", style = "margin-bottom: -10px;"),
            )
          })
          
          observeEvent(input$covariates,{
            if (input$covariates == "Covariate(s)") {
              
              shinyjs::show("covariates_variables")
              
              observeEvent(input$covariatesOptions,{
                if (!is.null(input$covariatesOptions)) {
                  output$chooseTest <- renderUI({
                    tagList(
                      
                      p(" ", style = "margin-top: 5px;"),
                      
                      h4(strong("Method?", style = "color:black")),
                      
                      p(" ", style = "margin-bottom: -20px;"),
                      
                      radioButtons("chooseMethod", label = "",#h4(strong("Method?", style = "color:black")),
                                   c("Linear regression"), selected = "Linear regression", width = '80%'),
                      div(id = "tauopt",
                          uiOutput("tauSlider")))
                  })
                } else {
                  output$chooseTest <- renderUI({
                    tagList(
                      #p(" ", style = "margin-bottom: -10px;"),
                      
                      p(" ", style = "margin-top: 5px;"),
                      
                      h4(strong("Method?", style = "color:black")),
                      
                      p(" ", style = "margin-bottom: -20px;"),
                      
                      radioButtons("chooseMethod", label = "",#h4(strong("Method?", style = "color:black")),
                                   c("Wilcoxon rank-sum test (default)", "Welch t-test"), selected = "Wilcoxon rank-sum test (default)", width = '80%'))
                  })
                }
              })
            } else if (input$covariates == "None") {
              shinyjs::hide("covariates_variables")
              output$chooseTest <- renderUI({
                tagList(
                  #p(" ", style = "margin-bottom: -10px;"),
                  
                  p(" ", style = "margin-top: 5px;"),
                  
                  h4(strong("Method?", style = "color:black")),
                  
                  p(" ", style = "margin-bottom: -20px;"),
                  
                  radioButtons("chooseMethod", label = "",#h4(strong("Method?", style = "color:black")), 
                               c("Wilcoxon rank-sum test (default)", "Welch t-test"), selected = "Wilcoxon rank-sum test (default)", width = '80%'))
              })
            }
          })
        }
      }
    })
    
    ################################################################################################
    
    ######################################
    ######################################
    ########## Beta diversity ############
    ######################################
    ######################################
    
    ######################################
    #### Cross-Sectional Data Analysis ###
    ######################################
    
    ## Beta Cross sectional UI options
    ## Beta primary variable UI option
    output$beta_primvar_cross <- renderUI({
      tagList(
        selectInput("beta.primvar_cross", label = h4(strong("Treatment Variable?", style = "color:black")),
                    c("Choose one" = "", sort(chooseData$nt.selected.bin)), selected = sort(chooseData$nt.selected.bin)[1], width = '80%'))
    })
    
    output$beta_prim_vars_types_cross <- renderUI({
      tagList(
        uiOutput("beta_morePrimvar_optcross"))
    })
    
    observeEvent(input$beta.primvar_cross,{
      
      if (input$beta.primvar_cross %in% chooseData$nt.selected.bin) {
        is.results$result = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.primvar_cross)
        
        nt.selected_cat =  chooseData$prim_varsSurv[[1]][ chooseData$prim_varsSurv[[1]] != input$beta.primvar_cross]
        chooseData$nt.selected.sub = nt.selected_cat
        
        if (is.results$result == "Binary") {
          beta.categos$cat1 = beta.bin.cat.func(chooseData$sam.dat, input$beta.primvar_cross)[1]
          beta.categos$cat2 = beta.bin.cat.func(chooseData$sam.dat, input$beta.primvar_cross)[2]
          
          output$beta_morePrimvar_optcross <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")),
              p("You can rename the categories of primary variable. MiSurv keeps up to 8 characters on graphs.", style = "font-size:11pt"), 
              textInput("betaCat1", label = (paste0("Reference: ",beta.categos$cat1)), value = beta.categos$cat1, width = '80%'),
              textInput("betaCat2", label = (paste0("Comparison: ",beta.categos$cat2)), value = beta.categos$cat2, width = '80%'))
          }) 
          
          #ntselected.prim_varsbin = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.primvar_cross)
          #ntselected.prim_vars = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.primvar_cross)
          nt.selected.prim <- c(chooseData$nt.selected.bin, chooseData$nt.selected.con)[!c(chooseData$nt.selected.bin, chooseData$nt.selected.con) %in% input$primvar.Surv.select]
          
          output$beta_covariates_cross <- renderUI({
            tagList(
              p(" ", style = "margin-bottom: -5px;"),
              
              prettyRadioButtons("beta_covariates_bin",label = h4(strong("Covariate(s)?", style = "color:black")), 
                                 status = "primary", icon = icon("check"), animation = "jelly",
                                 c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "beta_covariates_variables_bin", style = "margin-left: 2%",
                           prettyCheckboxGroup("beta_covariatesOptions_bin"," Please select covariate(s)", status = "primary",
                                               nt.selected.prim, width = '70%'))),
              
              #p(" ", style = "margin-bottom: -10px;"),
              
              p(" ", style = "margin-top: 5px;"),
              
              h4(strong("Method?", style = "color:black")),
              
              p(" ", style = "margin-bottom: -20px;"),
              
              radioButtons("beta_chooseMethod_bin", label = "",#h4(strong("Method?", style = "color:black")),
                           c("MiRKAT"), selected = "MiRKAT", width = '80%')#,
              
              # actionButton("beta_runbtn_cross_bin", (strong("Run!")), class = "btn-info")
            )
          })
          
          observeEvent(input$beta_covariates_bin,{
            if (input$beta_covariates_bin == "Covariate(s)") {
              shinyjs::show("beta_covariates_variables_bin")
            } else if (input$beta_covariates_bin == "None") {
              shinyjs::hide("beta_covariates_variables_bin")
            }
          })
        }
      }
    })
    
  })
  
  observeEvent(chooseData$taxa.out,{
    
    ################################################################################################
    
    ######################################
    ######################################
    ######## Taxonomic Analysis ##########
    ######################################
    ######################################
    
    output$primvars_taxa <- renderUI({
      tagList(
        prettyRadioButtons("dataType_taxa", label = h4(strong("Data Format?", style = "color:black")), animation = "jelly",
                           c("CLR (Default)", "Count (rarefied)", "Proportion", "Arcsine-root"), selected = "CLR (Default)",width = '70%'),
        p(" ", style = "margin-bottom: -10px;"),
        selectInput("primvar_taxa", label = h4(strong("Treatment Variable?", style = "color:black")),
                    choices = sort(chooseData$nt.selected.bin), selected = sort(chooseData$nt.selected.bin)[1], width = '80%'))
    })
    
    observeEvent(input$dataType_taxa,{
      if (input$dataType_taxa == "Count (rarefied)") {
        taxa.types$dataType = "rare.count"
        taxa.types$regression = "Negative binomial regression"
      } else if (input$dataType_taxa == "Proportion") {
        taxa.types$dataType = "imp.prop"
        taxa.types$regression = "Beta regression"
      } else if (input$dataType_taxa == "CLR (Default)") {
        taxa.types$dataType = "clr"
        taxa.types$regression = "Linear regression"
      } else if (input$dataType_taxa == "Arcsine-root") {
        taxa.types$dataType = "arcsin"
        taxa.types$regression = "Linear regression"
      }
    })
    
    observeEvent(input$primvar_taxa,{
      if (input$primvar_taxa %in% chooseData$prim_vars) {
        is.results$result = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar_taxa)
        
        nt.selected_cat =  chooseData$prim_varsSurv[[1]][ chooseData$prim_varsSurv[[1]] != input$primvar_taxa]
        chooseData$nt.selected.sub = nt.selected_cat
        
        if (is.results$result == "Binary") {
          taxa.categos$cat1 = taxa.bin.cat.func(chooseData$sam.dat, input$primvar_taxa)[1]
          taxa.categos$cat2 = taxa.bin.cat.func(chooseData$sam.dat, input$primvar_taxa)[2]
          
          output$morePrimvar_opt_taxa <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")), 
              p("You can rename the categories of primary variable. MiSurv keeps up to 8 characters on graphs.", style = "font-size:11pt"),
              textInput("taxaCat1", label = (paste0("Reference: ",taxa.categos$cat1)), value = taxa.categos$cat1, width = '80%'),
              textInput("taxaCat2", label = (paste0("Comparison: ",taxa.categos$cat2)), value = taxa.categos$cat2, width = '80%'))
          }) 
          
          nt.selected.prim <- c(chooseData$nt.selected.bin, chooseData$nt.selected.con)[!c(chooseData$nt.selected.bin, chooseData$nt.selected.con) %in% input$primvar.Surv.select]
          
          output$covariates_taxa <- renderUI({
            tagList(
              p(" ", style = "margin-bottom: -5px;"),
              prettyRadioButtons("covariates_taxa", label = h4(strong("Covariate(s)?", style = "color:black")), animation = "jelly",
                                 c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "covariates_variables_taxa", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariatesOptions_taxa"," Please select covariate(s)", status = "primary",
                                               nt.selected.prim, width = '70%'))),
              
              uiOutput("chooseTest_taxa"),
              
              p(" ", style = "margin-bottom: -10px;"),
              
              prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                 c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                 icon = icon("check"), width = '80%')
            )
          })
          
          observeEvent(input$covariates_taxa,{
            if (input$covariates_taxa == "Covariate(s)") {
              
              shinyjs::show("covariates_variables_taxa")
              
              observeEvent(input$covariatesOptions_taxa,{
                if (!is.null(input$covariatesOptions_taxa)) {
                  output$chooseTest_taxa <- renderUI({
                    tagList(
                      
                      p(" ", style = "margin-top: 5px;"),
                      
                      h4(strong("Method?", style = "color:black")),
                      
                      p(" ", style = "margin-bottom: -20px;"),
                      
                      radioButtons("chooseMethod_taxa", label = "", 
                                   c(taxa.types$regression), selected = taxa.types$regression,
                                   width = '80%'))
                  })
                } else {
                  
                  shinyjs::hide("covariates_variables_taxa")
                  output$chooseTest_taxa <- renderUI({
                    tagList(
                      p(" ", style = "margin-top: 5px;"),
                      
                      h4(strong("Method?", style = "color:black")),
                      
                      p(" ", style = "margin-bottom: -20px;"),
                      
                      radioButtons("chooseMethod_taxa", label = "", #h4(strong("Method?", style = "color:black")), 
                                   c("Wilcoxon rank-sum test (default)", "Welch t-test"), selected = "Wilcoxon rank-sum test (default)",)
                    )
                  })
                }
              })
            } else if (input$covariates_taxa == "None") {
              shinyjs::hide("covariates_variables_taxa")
              output$chooseTest_taxa <- renderUI({
                tagList(
                  
                  p(" ", style = "margin-top: 5px;"),
                  
                  h4(strong("Method?", style = "color:black")),
                  
                  p(" ", style = "margin-bottom: -20px;"),
                  
                  radioButtons("chooseMethod_taxa", label = "", #h4(strong("Method?", style = "color:black")), 
                               c("Wilcoxon rank-sum test (default)", "Welch t-test"),
                               selected = "Wilcoxon rank-sum test (default)", width = '80%'))
              })
            }
          }) 
        }
      }
    })
    
    ###############################################
    ######### Taxa Random Forest Analysis #########
    ###############################################
    
    output$model4_data_format <- renderUI({
      tagList(
        prettyRadioButtons("surv4.dataType_taxa", label = h4(strong("Data Format?", style = "color:black")), animation = "jelly",
                           c("CLR (Default)", "Count (rarefied)", "Proportion", "Arcsine-root"), selected = "CLR (Default)",width = '70%'))
    }) 
    
    observeEvent(input$surv4.dataType_taxa,{
      if (input$surv4.dataType_taxa == "Count (rarefied)") {
        taxa.types$dataType_model4 = "rare.count"           #module3 taxa count 수정 0812
      } else if (input$surv4.dataType_taxa == "Proportion") {
        taxa.types$dataType_model4 = "imp.prop"
      } else if (input$surv4.dataType_taxa == "CLR (Default)") {
        taxa.types$dataType_model4 = "clr"
      } else if (input$surv4.dataType_taxa == "Arcsine-root") {
        taxa.types$dataType_model4 = "arcsin"
      }
    })
    
    output$subgroupSurv.4 <- renderUI({
      tagList(
        
        h4(strong("Subgroup Analysis?", style = "color:black")), 
        p("You can select a subgroup of interest using a category of a nominal variable (usually the treatment variable) to perform subgroup analysis. If you select a subgroup, only the subjects in that subgroup are retained in the following analysis. Otherwise, all the subjects given in the data are retained (optional).", style = "font-size:11pt"),
        
        prettyRadioButtons("subgroup.4", label = NULL, status = "primary", icon = icon("check"),
                           animation = "jelly", choices = c("No", "Yes"), selected="No", width = '70%'),
        
        p(" ", style = "margin-bottom: -25px;"),
        
        shinyjs::hidden(
          tagList(
            shiny::div(id = "subgroup.sel.4", style = "margin-left: 5%",
                       selectInput("subgroup.sel.4", label = "", 
                                   c("Choose one" = "", "null"), selected = "null", width = '70%'))
          )
        ),
        p(" ", style = "margin-bottom: -10px;"),
        uiOutput("pick_subgroup.4")
      )
    })
    
    
    observeEvent(input$subgroup.4,{
      if (input$subgroup.4 == "Yes"){
        cat.for.subgroup = chooseData$prim_varsCox[[1]][ !chooseData$prim_varsCox[[1]]%in% input$censor.select ]
        
        if(!is.null(input$subgroup.sel.4)){
          if(input$subgroup.sel.4 %in% cat.for.subgroup) {
            updateSelectInput(session, "subgroup.sel.4", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                              selected = input$subgroup.sel.4)
          } else {
            updateSelectInput(session, "subgroup.sel.4", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                              selected = sort(cat.for.subgroup)[1])
          } 
        } else {
          updateSelectInput(session, "subgroup.sel.4", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                            selected = sort(cat.for.subgroup)[1])
        }
        shinyjs::show("subgroup.sel.4")
        shinyjs::show("pick_subgroup.4")
      }
      else {
        shinyjs::hide("subgroup.sel.4")
        shinyjs::hide("pick_subgroup.4")
      }
    })
    
    
    observeEvent(input$subgroup.sel.4, {
      if(input$subgroup.sel.4 %in% infile$ori.var[[1]]){
        
        subgroup.val <- alpha.bin.cat.ref.ori.func(chooseData$sam.dat, input$subgroup.sel.4)
        
        output$pick_subgroup.4 <- renderUI({
          tagList(
            p(" ", style = "margin-bottom: -5px;"),
            shiny::div(id = "pick_subgroup.4", style = "margin-left: 5%",
                       radioButtons("pick_subgroup.4", label =  "Please select a subgroup:", 
                                    choices = subgroup.val,
                                    selected = subgroup.val[1], width = '70%')),
            p(" ", style = "margin-bottom: -20px;"),
          )
        })
      }
    })
    
    output$method.4 <- renderUI({
      tagList(
        p(" ", style = "margin-top: +40px;"),
        
        h4(strong("Method?", style = "color:black")),
        
        p(" ", style = "margin-bottom: -20px;"),
        
        radioButtons("surv4.method.select", label = "",
                     c("Random survival forests"), selected = "Random survival forests", width = '70%'),
        
        p(" ", style = "margin-bottom: -10px;"),
        
        shinyjs::hidden(
          shiny::div(id = "surv4.fold.select", #style = "margin-left: 2%",
                     selectInput("surv4.fold.select", label = h4(strong("Number of Folds?", style = "color:black")),
                                 c("Choose one" = "", c(5, 10)), selected = 10, width = '70%'))),
        uiOutput("surv4.trees"),
        uiOutput("surv4.num.display"),
        
        p(" ", style = "margin-bottom: -10px;"),
        
        prettyRadioButtons("include_speciesM4", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                           c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                           icon = icon("check"), width = '80%')
        
      )
    })
    
    observeEvent(input$surv4.method.select, {
      if (input$surv4.method.select == "Random survival forests") {
        shinyjs::show("surv4.trees")
        shinyjs::hide("surv4.fold.select")
        output$surv4.trees <- renderUI({
          tagList(
            p(" ", style = "margin-top: +30px;"),
            h4(strong("# of Trees?", style = "color:black")),
            p("The number of trees to be ensembled (default: 5,000).", style = "font-size:11pt"),
            selectInput("surv4.tree.select", label = NULL,
                        c("Choose one" = "", c(5000, 10000)), selected = "5000", width = '70%')
          )
        })
        shinyjs::show("surv4.num.display")
        output$surv4.num.display <- renderUI({
          tagList(
            p(" ", style = "margin-top: 5px;"),
            p(strong("# of Taxa?"), style = "font-size:13.5pt; color:black"),
            p(" ", style = "margin-bottom: -5px;"),
            Num_PredComment,
            p(" ", style = "margin-bottom: -20px;"),
            selectInput("surv4.num.display", label = "",
                        c("Choose one" = "", c(10, 20, 30)), selected = 10, width = '70%')
          )
        })
      } else {
        shinyjs::hide("surv4.trees")
        shinyjs::hide("surv4.num.display")
        shinyjs::show("surv4.fold.select")
      }
    })
    
    
    ######################################
    ######### Taxa Surv Analysis #########
    ######################################
    
    output$censor.selectCoxT <- renderUI({
      tagList(
        prettyRadioButtons("surv.dataType_taxa", label = h4(strong("Data Format?", style = "color:black")), animation = "jelly",
                           c("CLR (Default)", "Count (rarefied)", "Proportion", "Arcsine-root"), selected = "CLR (Default)",width = '70%'),
        uiOutput("surv3T.subgroup")
      )
    })
    
    output$subgroupSurv.T <- renderUI({
      tagList(
        
        h4(strong("Subgroup Analysis?", style = "color:black")), 
        p("You can select a subgroup of interest using a category of a nominal variable (usually the treatment variable) to perform subgroup analysis. If you select a subgroup, only the subjects in that subgroup are retained in the following analysis. Otherwise, all the subjects given in the data are retained (optional).", style = "font-size:11pt"),
        
        prettyRadioButtons("subgroup.t", label = NULL, status = "primary", icon = icon("check"),
                           animation = "jelly", choices = c("No", "Yes"), selected="No", width = '70%'),
        p(" ", style = "margin-bottom: -25px;"),
        shinyjs::hidden(
          tagList(
            shiny::div(id = "subgroup.sel.t", style = "margin-left: 5%",
                       selectInput("subgroup.sel.t", label = "", 
                                   c("Choose one" = "", "null"), selected = "null", width = '70%'))
          )
        ),
        p(" ", style = "margin-bottom: -10px;"),
        uiOutput("pick_subgroup.t"),
        p(" ", style = "margin-bottom: +20px;"),
      )
    })
    
    observeEvent(input$subgroup.t,{
      if (input$subgroup.t == "Yes"){
        cat.for.subgroup = chooseData$prim_varsCox[[1]][ !chooseData$prim_varsCox[[1]]%in% input$censor.select ]
        
        if(!is.null(input$subgroup.sel.t)){
          if(input$subgroup.sel.t %in% cat.for.subgroup) {
            updateSelectInput(session, "subgroup.sel.t", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                              selected = input$subgroup.sel.t)
          } else {
            updateSelectInput(session, "subgroup.sel.t", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                              selected = sort(cat.for.subgroup)[1])
          } 
        } else {
          updateSelectInput(session, "subgroup.sel.t", label = "", c("Choose one" = "", sort(cat.for.subgroup)), 
                            selected = sort(cat.for.subgroup)[1])
        }
        shinyjs::show("subgroup.sel.t")
        shinyjs::show("pick_subgroup.t")
      }
      else {
        shinyjs::hide("subgroup.sel.t")
        shinyjs::hide("pick_subgroup.t")
      }
    })
    
    observeEvent(input$subgroup.sel.t, {
      if(input$subgroup.sel.t %in% infile$ori.var[[1]]){
        
        subgroup.val <- alpha.bin.cat.ref.ori.func(chooseData$sam.dat, input$subgroup.sel.t)
        
        output$pick_subgroup.t <- renderUI({
          tagList(
            p(" ", style = "margin-bottom: -5px;"),
            shiny::div(id = "pick_subgroup.t", style = "margin-left: 5%",
                       radioButtons("pick_subgroup.t", label =  "Please select a subgroup:", 
                                    choices = subgroup.val,
                                    selected = subgroup.val[1], width = '70%')),
            p(" ", style = "margin-bottom: -20px;")
          )
        })
      }
    })
    
    observeEvent(input$surv.dataType_taxa,{
      if (input$surv.dataType_taxa == "Count (rarefied)") {
        taxa.types$dataType_Surv = "rare.count"           #module3 taxa count 수정 0812
      } else if (input$surv.dataType_taxa == "Proportion") {
        taxa.types$dataType_Surv = "imp.prop"
      } else if (input$surv.dataType_taxa == "CLR (Default)") {
        taxa.types$dataType_Surv = "clr"
      } else if (input$surv.dataType_taxa == "Arcsine-root") {
        taxa.types$dataType_Surv = "arcsin"
      }
    })
    
    observeEvent(input$covariatesCoxT,{
      if (input$covariatesCoxT == "Covariate(s)") {
        
        shinyjs::show("covariates_variablesCoxT")
        
      } else if (input$covariatesCoxT == "None") {
        
        shinyjs::hide("covariates_variablesCoxT")
        
      }
    })
    
    output$covariatesCoxT <- renderUI({
      tagList(
        p(" ", style = "margin-top: +25px;"),
        
        prettyRadioButtons("covariatesCoxT",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                           animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
        
        shinyjs::hidden(
          shiny::div(id = "covariates_variablesCoxT", style = "margin-left: 2%",
                     prettyCheckboxGroup("covariatesOptionsCoxT"," Please select covariate(s)", status = "primary",
                                         c(chooseData$nt.selected.bin, chooseData$nt.selected.con), width = '70%'))),
        
        #p(" ", style = "margin-bottom: -10px;"),
        
        p(" ", style = "margin-top: 5px;"),
        
        h4(strong("Method?", style = "color:black")),
        
        p(" ", style = "margin-bottom: -20px;"),
        
        radioButtons("surv3.taxa.method", label = "",#h4(strong("Method?", style = "color:black")),
                     c("Cox model"), selected = "Cox model", width = '80%'),
        
        p(" ", style = "margin-bottom: -10px;"),
        
        prettyRadioButtons("include_species_model3", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                           c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                           icon = icon("check"), width = '80%')#,
        # actionButton("runbtn_CoxT", (strong("Run!")), class = "btn-info")
      )
    })
  })
  
  #########################################################################################################
  
  ######################################
  ######################################
  ##########  RUN BUTTONS   ############
  ######################################
  ######################################
  
  ##########################################
  ##     Survival Data Analysis           ##
  ##########################################
  
  observeEvent(input$surv.KMrun, {
    validate(
      if (input$covariates.surv == "Covariate(s)" & is.null(input$covariatesOptions.surv)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("surv.KMrun")
    # shinyjs::disable("chooseAdjustment")
    # shinyjs::disable("primvar")
    # shinyjs::disable("chooseMethod")
    # shinyjs::disable("covariates")
    # shinyjs::disable("covariatesOptions")
    # shinyjs::disable("alphaCat1")
    # shinyjs::disable("alphaCat2")
    # shinyjs::disable("chooseData")
    
    withProgress({
      
      #########################################################
      incProgress(1/10, message = "Applying changes")
      
      rename.cats_ref <- input$prim.var.rename1.Surv
      rename.cats_com <- input$prim.var.rename2.Surv
      
      Surv.bin.cat.ref.ori.out <- alpha.bin.cat.ref.ori.func(chooseData$sam.dat, input$primvar.Surv.select)
      
      sam_dat0 <- alpha.bin.cat.recode.func(chooseData$sam.dat, input$primvar.Surv.select, Surv.bin.cat.ref.ori.out,
                                            rename.cats_ref, rename.cats_com)
      
      print("here1")
      
      sample_here <<- sam_dat0
      covariates_here <<- input$covariatesOptions.surv
      
      sam_dat0 <- sam_no_miss_cov(sam_dat0, input$covariatesOptions.surv)
      print("here2")
      print(sam_dat0)
      
      # sam_dat0 <- sam.follow.time.func(sam_dat, chooseData$surv.Time.select, chooseData$censor.select, chooseData$follow.time)
      
      # if (chooseData$subgroup == "Yes"){
      #   sam_dat0 <- sam_dat0[ sam_dat0[[chooseData$subgroup.sel]] == chooseData$pick_subgroup ]
      # }
      
      #########################################################
      incProgress(3/10, message = "Displaying Variables in Use")
      
      #num.treat1_size <- surv.num.treat1(sam_dat, input$primvar.Surv.select)
      #num.treat2_size <- surv.num.treat2(sam_dat, input$primvar.Surv.select)
      
      surv.dat0 <- surv.events.dat.func(sam_dat0, 
                                        input$primvar.Surv.select, 
                                        chooseData$surv.Time.select, 
                                        chooseData$censor.select)
      
      num.treat1_size <- surv.num.treat1(surv.dat0, "treat")
      num.treat2_size <- surv.num.treat2(surv.dat0, "treat")
      
      num.event <- surv.num.uncen(surv.dat0)
      num.cen <- surv.num.censr(surv.dat0)
      
      output$treatment_size <- renderValueBox({
        treat.ls <- unique( sam_dat0[[input$primvar.Surv.select]] )
        
        valueBox(
          value = tags$p(paste0(num.treat1_size), style = "font-size: 75%;"),
          paste("Sample Size (", treat.ls[1],")", sep="" ), icon = icon("user-circle"), color = "red")
      })
      
      output$control_size <- renderValueBox({
        treat.ls <- unique( sam_dat0[[input$primvar.Surv.select]] )
        
        valueBox(
          value = tags$p(paste0(num.treat2_size), style = "font-size: 75%;"),
          paste("Sample Size (", treat.ls[2],")", sep="" ), icon = icon("user-circle"), color = "orange")
      })
      
      output$censored_size <- renderValueBox({
        
        valueBox(
          value = tags$p(paste0(num.cen), style = "font-size: 75%;"),
          "Sample Size (Censored)", icon = icon("user-circle"), color = "yellow")
      })
      
      output$uncensored_size <- renderValueBox({
        
        valueBox(
          value = tags$p(paste0(num.event), style = "font-size: 75%;"),
          "Sample Size (Event)", icon = icon("user-circle"), color = "teal")
      })
      
      
      validate(
        if (num.treat1_size < 2 & num.treat2_size < 2 & num.event < 2 ) {
          if( num.treat1_size < 2 & num.treat2_size < 2 ){
            showNotification("Error: Too few sample Treatment size. Select different Treatment variable.",
                             type = "error")
          }
          if( num.event < 2 ){
            showNotification("Error: Too few sample Event size. Select different Event variable.",
                             type = "error")
          }
        }
      )
      
      #########################################################
      incProgress(6/10, message = "Calculating")
      
      # If Covariates = None, usual ( but if Cox is checked, get cox curve)
      # If Covariates!= None, get Cox PH
      
      #sam_dat <- resue.sam.dat #chooseData$sam_dat
      
      if ( input$covariates.surv == "None") {
        
        surv.dat <- surv.events.dat.func(sam_dat0, 
                                         input$primvar.Surv.select, 
                                         chooseData$surv.Time.select, 
                                         chooseData$censor.select)
        
        if ( input$pick_KM_Cox == "Cox model") {
          cox.fit <- surv.Cox.fit.func1( surv.dat )
          
          output$surv_display_results = renderUI({
            tagList(
              box(title = strong( "Survival Curve"), 
                  align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                  plotOutput("surv.plot2", height = 400, width = 500)
              )
            )
          })
          output$surv.plot2 = renderPlot({
            isolate(surv.Cox.plot.func( cox.fit, surv.dat ))
          })
        }
        else {
          surv.KM.fit <- surv.KM.fit.func( surv.dat )
          
          output$surv_display_results = renderUI({
            tagList(
              box(title = strong( "Kaplan-Meier Curve"), 
                  align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                  plotOutput("surv.plot3", height = 400, width = 500)
              )
            )
          })
          output$surv.plot3 = renderPlot({
            isolate(surv.KM.plot.func( surv.KM.fit, surv.dat, input$sig.test.Surv ))
          })
        }
      }
      else {
        
        surv.dat <- surv.events.dat.func2(sam_dat0, 
                                          input$primvar.Surv.select, 
                                          chooseData$surv.Time.select, 
                                          chooseData$censor.select,
                                          input$covariatesOptions.surv)
        
        cox.model <- surv.Cox.fit.func2( surv.dat, input$covariatesOptions.surv)
        
        output$surv_display_results = renderUI({
          tagList(
            box(title = strong( "Adjusted Survival Curve"),
                align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("surv.plot1", height = 400, width = 500)))
        })
        
        output$surv.plot1 = renderPlot({
          isolate(surv.Cox.plot.func( cox.model, surv.dat, input$covariatesOptions.surv))
        })
      }
      
      if( input$pick_KM_Cox == "Cox model"){
        ref_string = REFERENCE_CHECK(method_name = isolate(input$pick_KM_Cox) )
        
      }
      else{
        ref_string = REFERENCE_CHECK(method_name = isolate(input$pick_KM_Cox), sig_test = isolate(input$sig.test.Surv) )
      }
      
      if (is.null(ref_string)) {
        shinyjs::hide("referencesM1")
      } else {
        shinyjs::show("referencesM1")
        output$referencesM1 = renderUI({
          tagList(
            box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                HTML(paste(ref_string, collapse="<br/>"))
            )
          )
        })
      }
      shinyjs::enable("surv.KMrun")
      # shinyjs::enable("primvar")
      # shinyjs::enable("chooseAdjustment")
      # shinyjs::enable("chooseMethod")
      # shinyjs::enable("covariates")
      # shinyjs::enable("covariatesOptions")
      # shinyjs::enable("alphaCat1")
      # shinyjs::enable("alphaCat2")
      # shinyjs::enable("chooseData")
    })
    
  })
  
  ##########################################
  ## Alpha Cox Proportion Hazard Analysis ##
  ##########################################
  observeEvent(input$runbtn_CoxA, {
    validate(
      if (input$covariatesCoxA == "Covariate(s)" & is.null(input$covariatesOptionsCoxA)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtn_CoxA")
    # shinyjs::disable("chooseAdjustment")
    # shinyjs::disable("primvar")
    # shinyjs::disable("chooseMethod")
    # shinyjs::disable("covariates")
    # shinyjs::disable("covariatesOptions")
    # shinyjs::disable("alphaCat1")
    # shinyjs::disable("alphaCat2")
    # shinyjs::disable("chooseData")
    
    withProgress(
      
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        ##############################################
        incProgress(1/10, message = "Calculating")
        
        #chooseData$alpha.div = chooseData$alpha.div.rare
        
        ##############################################
        incProgress(5/10, message = "Plotting")
        
        new_data <- data.union.func(chooseData$sam.dat, chooseData$alpha.div)
        
        new_sam.dat   <- new_data$new_sam.dat
        new_alpha.div <- new_data$new_alpha.div
        
        
        if (input$subgroup.a== "Yes"){
          ind <- new_sam.dat [[input$subgroup.sel.a]] == input$pick_subgroup.a
          new_sam.dat  <- new_sam.dat [ind, ]
          new_alpha.div <- new_alpha.div [ind,]
        }
        
        new_alpha.div <- alpha_no_miss_cov(new_sam.dat,  new_alpha.div, input$covariatesOptionsCoxA)
        new_sam.dat <- sam_no_miss_cov(new_sam.dat, input$covariatesOptionsCoxA)
        
        new_surv.dat <- alpha.surv.events.dat.func(new_sam.dat, chooseData$surv.Time.select, chooseData$censor.select)
        
        if (input$covariatesCoxA =="None") {
          alpha.cox.out <- alpha.cox.test(new_surv.dat, new_alpha.div)
        }
        else {
          cov.dat <- subset(new_sam.dat, select = input$covariatesOptionsCoxA)
          alpha.cox.out <- alpha.cox.test(new_surv.dat, new_alpha.div, cov.dat)
        }
        
        multi.test <- FALSE
        alphaS.down.results$alpha.cox.out = alpha.cox.out
        
        # if (input$chooseAdjustment3A == "Yes") {
        #   multi.test <- TRUE
        #   alphaS.down.results$alpha.cox.out = q.func(alpha.cox.out, method = "BH")
        # }
        # else if (input$chooseAdjustment3A == "No (Default)") {
        #   multi.test <- FALSE
        #   alphaS.down.results$alpha.cox.out = alpha.cox.out
        # }
        
        alphaS.down.results$alpha.cox.forest.plot <- alpha.cox.forest.plot(alphaS.down.results$alpha.cox.out, multi.test = multi.test)
        
        output$alpha_surv_display_resultsCoxA = renderUI({
          tagList(
            box(title = strong( "Graphs for Alpha Diversity"),
                align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("surv.plot.coxA", height = 400, width = 500)))
        })
        
        output$surv.plot.coxA = renderPlot({
          alphaS.down.results$alpha.cox.forest.plot
        })
        
        output$alpha_surv_downloadTable = renderUI({
          tagList(
            p(" ", style = "margin-top: 20px;"),
            box(title = strong("Download Output Table"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Data Analysis Outputs"),
                downloadButton("downloadTabl_aS1", "Download", width = '50%', style = "background-color: red3"), br(),
                h5("Data Analysis Plot Image"),
                downloadButton("downloadTabl_aS2", "Download", width = '50%', style = "background-color: red3")
            )
          )
        })
        output$downloadTabl_aS1 <- downloadHandler(
          filename = function() {
            paste("Alpha.Survival.Table.csv")
          },
          content = function(file) {
            write.csv(alphaS.down.results$alpha.cox.out, file)
          }
        )
        output$downloadTabl_aS2 <- downloadHandler(
          
          filename = function() {
            paste("Alpha.Survival.Output.png")
          },
          content = function(file) {
            #write.csv(alphaS.down.results$alpha.cox.forest.plot, file)
            png(file=file,
                width=650, height=400)
            print( alphaS.down.results$alpha.cox.forest.plot )
            dev.off()
          }
        )
        
        ref_string = REFERENCE_CHECK(method_name = isolate(input$surv3.alpha.method), FDR = "No (Default)")
        
        if (is.null(ref_string)) {
          shinyjs::hide("referencesM3.alpha")
        } else {
          shinyjs::show("referencesM3.alpha")
          output$referencesM3.alpha = renderUI({
            tagList(
              box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
      }
    )
    shinyjs::enable("runbtn_CoxA")
  })
  
  ##########################################
  ## Beta Cox Proportion Hazard Analysis ##
  ##########################################
  observeEvent(input$runbtn_CoxB, {
    print("hello come")
    validate(
      if (input$covariatesCoxB == "Covariate(s)" & is.null(input$covariatesOptionsCoxB)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    shinyjs::disable("runbtn_CoxB")
    
    
    withProgress(
      
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        ##############################################
        
        incProgress(1/10, message = "Calculating")
        
        betaS.primvar_time = chooseData$surv.Time.select
        betaS.primvar_cross  = chooseData$censor.select
        
        beta.bin.cat.ref.ori.out <- beta.bin.cat.ref.ori.func(chooseData$sam.dat, betaS.primvar_cross)
        
        rename.catsbin_ref = beta.bin.cat.ref.ori.out[1]
        rename.catsbin_com = beta.bin.cat.ref.ori.out[2]
        
        beta.sam_dat.bin <- beta.bin.cat.recode.func(chooseData$sam.dat, betaS.primvar_cross,
                                                     beta.bin.cat.ref.ori.out,
                                                     rename.catsbin_ref, rename.catsbin_com)
        
        beta.sam_dat.bin <- sam_no_miss_cov(beta.sam_dat.bin, input$covariatesOptionsCoxB)
        
        ds.Ks.res <- ds.Ks$res
        ds.Ks.res <- beta_no_miss_cov(beta.sam_dat.bin, ds.Ks.res, input$covariatesOptionsCoxB)
        
        if (input$subgroup.b == "Yes"){
          ind.sub <- beta.sam_dat.bin [[input$subgroup.sel.b]] == input$pick_subgroup.b
          beta.sam_dat.bin  <- beta.sam_dat.bin [ind.sub, ]
        }
        
        
        for(i in 1:2){
          for(j in 1:length(ds.Ks$res[[i]])) {
            ind <- rownames(ds.Ks$res[[i]][[j]]) %in% rownames(beta.sam_dat.bin)
            ds.Ks.res[[i]][[j]] <- ds.Ks$res[[i]][[j]][ind,ind]
          }
        }
        
        if (input$covariatesCoxB == "None") {
          betaS.bin.out <- betaS.bin.cat.ref.func(betaS.primvar_cross, betaS.primvar_time,
                                                  rename.catsbin_ref, rename.catsbin_com,
                                                  beta.sam_dat.bin, Ds.Ks = ds.Ks.res)
          betaS.data.results$data.q.out <- betaS.bin.out
        } 
        
        else if (input$covariatesCoxB == "Covariate(s)") {
          if (is.null(input$covariatesOptionsCoxB)) {
            betaS.bin.out <- betaS.bin.cat.ref.func(betaS.primvar_cross, betaS.primvar_time,
                                                    rename.catsbin_ref, rename.catsbin_com,
                                                    beta.sam_dat.bin, Ds.Ks = ds.Ks.res)
            betaS.data.results$data.q.out <- betaS.bin.out
          } 
          else {
            betaS.bin.cov.out <- betaS.bin.cov.cat.ref.func(betaS.primvar_cross, betaS.primvar_time,
                                                            rename.catsbin_ref, rename.catsbin_com,
                                                            input$covariatesOptionsCoxB, beta.sam_dat.bin,
                                                            Ds.Ks = ds.Ks.res)
            betaS.data.results$data.q.out <- betaS.bin.cov.out
          }
        }
        
        if (isolate(input$covariatesCoxB) == "None") {
          
          betaS.down.results$CS = mirkatS.bin(betaS.data.results$data.q.out)
          betaS.plot.info <- mirkatS.bin.plot1(betaS.down.results$CS, betaS.data.results$data.q.out)
          
          output$m3_beta_2d = renderPlot({
            isolate(try(mirkatS.bin.plot2(betaS.down.results$CS, betaS.data.results$data.q.out, betaS.plot.info$mod, betaS.plot.info$sub.tit), silent = TRUE))
          })
        } 
        
        else if (isolate(input$covariatesCoxB) == "Covariate(s)") {
          
          betaS.down.results$CS = mirkatS.bin.cov(betaS.data.results$data.q.out)
          
          betaS.plot.info <- mirkatS.bin.plot1(betaS.down.results$CS, betaS.data.results$data.q.out)
          
          output$m3_beta_2d = renderPlot({
            isolate(try( mirkatS.bin.plot2(betaS.down.results$CS, betaS.data.results$data.q.out, betaS.plot.info$mod, betaS.plot.info$sub.tit), silent = TRUE))
          })
        }
        
        betaS.data.results$betaS.plot.info <- isolate(betaS.plot.info)
        
        # output$betaS_display_results_cross = renderUI({
        #   tagList(
        #     tabPanel(title = strong("Graphs for Beta Diversity"), align = "center",
        #              tabsetPanel(id = "m3betaTabs",
        #                          tabPanel(title = "2D PCoA", align = "center",
        #                                   plotOutput("m3_beta_2d", height = 850, width = 500)))
        #     )
        #   )
        # })
        
        output$betaS_display_results_cross = renderUI({
          tagList(
            box(title = strong("Graphs for Beta Diversity"), 
                align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("m3_beta_2d", height = 850, width = 500)
            )
          )
        })
        
        beta_ind <- c("Jaccard", "Bray-Curtis", "U.UniFrac", "G.UniFrac", "W.UniFrac")
        
        output$beta_surv_downloadTable = renderUI({
          tagList(
            p(" ", style = "margin-top: 20px;"),
            box(title = strong("Download Output Table"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Data Analysis Outputs"),
                downloadButton("downloadTabl_bS1", "Download", width = '50%', style = "background-color: red3"), br(),
                h5("Data Analysis Plot Image"),
                downloadButton("downloadTabl_bS2", "Download", width = '50%', style = "background-color: red3")
            )
          )
        })
        
        output$downloadTabl_bS1 <- downloadHandler(
          filename = function() {
            paste("Beta.Survival.Table.csv")
          },
          content = function(file) {
            out_temp = as.data.frame(unlist(betaS.down.results$CS))
            rownames(out_temp) = c("Jaccard","Bray.Curtis","U.UniFrac","G.UniFrac","W.UniFrac", "OMiRKAT-S_p")
            colnames(out_temp) = "p-value"
            betaS.down.results$CS = out_temp
            #write.table(betaS.down.results$CS, file, sep="\t")
            write.csv(betaS.down.results$CS, file)
          }
        )
        output$downloadTabl_bS2 <- downloadHandler(
          filename = function() {
            paste("Beta.Survival.Output.png")
          },
          content = function(file) {
            png(file=file,
                width=1050, height=950)
            print( mirkatS.bin.plot2(betaS.down.results$CS, betaS.data.results$data.q.out, betaS.plot.info$mod, betaS.plot.info$sub.tit) )
            dev.off()
          }
        )
        ref_string = REFERENCE_CHECK(method_name = isolate(input$surv3.beta.method), FDR = "No")
        if (is.null(ref_string)) {
          shinyjs::hide("referencesM3.beta")
        } 
        else {
          shinyjs::show("referencesM3.beta")
          output$referencesM3.beta = renderUI({
            tagList(
              box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
      }
    )
    
    shinyjs::enable("runbtn_CoxB")
    
  })
  ##########################################
  ## Taxa Cox Proportion Hazard Analysis ###
  ##########################################
  observeEvent(input$runbtn_CoxT, {
    validate(
      if (input$covariatesCoxT == "Covariate(s)" & is.null(input$covariatesOptionsCoxT)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtn_CoxT")
    
    
    
    taxa.sam.dat <- chooseData$sam.dat[match(rownames(chooseData$taxa.out$clr$phylum),rownames(chooseData$sam.dat)),]
    taxa.out.surv <- chooseData$taxa.out[[taxa.types$dataType_Surv]]
    taxa.out.surv <- taxa_no_miss_cov_ind(taxa.sam.dat, taxa.out.surv, input$covariatesOptionsCoxT)
    
    taxa.names.out.surv <- chooseData$taxa.names.out
    
    if (input$subgroup.t== "Yes"){
      ind.sub <- taxa.sam.dat[[input$subgroup.sel.t]] == input$pick_subgroup.t
      taxa.sam.dat <- taxa.sam.dat[-ind.sub,]  
      taxa.out.surv <- lapply(taxa.out.surv, function(x){
        return(x[-ind.sub,])})
    }
    
    validate(
      if (sum(taxa.sam.dat[,chooseData$censor.select], na.rm = T) <= 0) {
        showNotification("Error: No Events in Data",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtn_CoxT")
    
    withProgress(
      
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        ##############################################
        incProgress(1/10, message = "Calculating")
        
        #chooseData$alpha.div = chooseData$alpha.div.rare
        
        ##############################################
        incProgress(5/10, message = "Plotting")
        
        if (input$include_species_model3 == "Phylum - Species") {
          include = TRUE
        } else {
          include = FALSE
        }
        
        incProgress(1/10, message = "Examining Data in progress")
        
        survtime <- taxa.sam.dat[,chooseData$surv.Time.select]
        status <- taxa.sam.dat[,chooseData$censor.select]
        
        if (input$covariatesCoxT =="None") {
          taxa.cox.out <- taxa.cox.test(survtime, status, taxa.out.surv)
        }
        else {
          taxa.cox.out <- taxa.cox.test(survtime, status, taxa.out.surv, cov.dat = taxa.sam.dat[,input$covariatesOptionsCoxT])
        }
        
        nrow <- surv.taxa.forest.plot.pages(taxa.cox.out, include, mult.test.cor = "TRUE")
        forestplot.data <- surv.taxa.forest.plot.pages1(taxa.cox.out, taxa.names.out.surv, report.type = "HR", species.include = include, mult.test.cor = TRUE)
        
        if (any(!is.na(unlist(taxa.names.out.surv$duplicates)))) {
          duplicate.taxa <- sapply(strsplit(unlist(taxa.names.out.surv$duplicates), " :"),  "[", 1)
          taxon.inplot <- unlist(lapply(forestplot.data$all.text.tab, `[`, i =, j = 3))
          duplicate.texts <- sum(duplicate.taxa %in% taxon.inplot)
        } else {
          duplicate.texts <- 0
        }
        
        if (duplicate.texts > 0) {
          output$surv_taxa_results = renderUI({
            tagList(
              do.call(tabsetPanel, lapply(1:nrow, function(i) {
                tabPanel(title = paste0("Page ", i), align = "center",
                         plotOutput(paste0("forest_Surv", i), height = 800, width = 750),
                         plotOutput(paste0("duplicates_surv", i), height = 11.5*duplicate.texts+10, width = 750))
              })) 
            )
          })
        } else {
          output$surv_taxa_results= renderUI({
            tagList(
              do.call(tabsetPanel, lapply(1:nrow, function(i) {
                tabPanel(title = paste0("Page ", i), align = "center",
                         plotOutput(paste0("forest_Surv", i), height = 800, width = 750))
              })) 
            )
          })
        }
        
        lapply(1:nrow, function(j) {
          output[[paste0("forest_Surv", j)]] <- renderPlot({
            taxa.surv.forest.plot.pages2(page.taxa.q.out = forestplot.data, page = j)
          })
        })
        
        lapply(1:nrow, function(j) {
          output[[paste0("duplicates_surv", j)]] <- renderPlot({
            duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
          })
        })
        
        taxa.sig <- taxa.sig.dend(taxa.cox.out, chooseData$NAadded$tax.tab, "twopi", include)
        
        flow.text <- taxa.sig$flow.text
        taxon.tab <- taxa.sig$taxon.tab
        ci.tab.all <- taxa.sig$ci.tab.all
        
        if ( include == FALSE){
          taxon.tab <- taxon.tab[ !grepl("s_", taxon.tab$Taxon, fixed = TRUE), ]
        }
        
        if ( length(ci.tab.all) > 1 ){
          
          for( i in 1:nrow(taxon.tab)){
            if ( ci.tab.all[-1][i] < 0){
              taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "blue")
              taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
            }
            else{
              taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "red")
              taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
            }
          }
        }
        
        N <- dim(taxon.tab)[1]
        itr <- ceiling(N/5)
        tab.one   <- data.frame( matrix(ncol=2,nrow=0) )
        tab.two   <- data.frame( matrix(ncol=2,nrow=0) )
        tab.three <- data.frame( matrix(ncol=2,nrow=0) )
        tab.four  <- data.frame( matrix(ncol=2,nrow=0) )
        tab.five  <- data.frame( matrix(ncol=2,nrow=0) )
        
        colnames(tab.one) <- c("ID", "Taxon")
        colnames(tab.two) <- c("ID", "Taxon")
        colnames(tab.three) <- c("ID", "Taxon")
        colnames(tab.four) <- c("ID", "Taxon")
        colnames(tab.five) <- c("ID", "Taxon")
        
        if ( dim(taxon.tab)[1] > 0 ) {
          
          i = 1
          j = 1
          N.rmd <- N %% 5
          N.fix <- N + 5 - N.rmd
          while ( i <= itr) {
            tab.one  [i,] <- taxon.tab[j,]
            tab.two  [i,] <- taxon.tab[j+1,]
            tab.three[i,] <- taxon.tab[j+2,]
            tab.four [i,] <- taxon.tab[j+3,]
            tab.five [i,] <- taxon.tab[j+4,]
            i <- i + 1  
            j <- j + 5
          }
          row.names(tab.one) <- NULL
          row.names(tab.two) <- NULL
          row.names(tab.three) <- NULL
          row.names(tab.four) <- NULL
          row.names(tab.five) <- NULL
          
          tab.one <- na.omit(tab.one)
          tab.two <- na.omit(tab.two)
          tab.three <- na.omit(tab.three)
          tab.four <- na.omit(tab.four)
          tab.five <- na.omit(tab.five)
        }
        
        output$survival_taxa_display_dend = renderUI({
          
          box(title = strong("Dendrogram"), width = 12, status = "primary", solidHeader = TRUE,
              
              fluidRow(width = 12, align = "center",
                       div(style = "display: inline-block:vertical-align:top;", grVizOutput("dendrogram_surv", height = 1000, width = 1000)) ),
              br(),
              fluidRow(width = 12, align = "center",
                       tagList(
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("M3sig_1_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("M3sig_2_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("M3sig_3_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("M3sig_4_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("M3sig_5_taxonlist") )
                       )
              )
          )
        })
        
        output$dendrogram_surv = renderGrViz({
          flow.text
        })
        output$M3sig_1_taxonlist <- renderText({
          sig.tab1 <- kable(tab.one, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab1
        })
        output$M3sig_2_taxonlist <- renderText({
          
          sig.tab2 <- kable(tab.two, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab2
        })
        output$M3sig_3_taxonlist <- renderText({
          
          sig.tab3 <- kable(tab.three, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab3
        })
        output$M3sig_4_taxonlist <- renderText({
          
          sig.tab4 <- kable(tab.four, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab4
        })
        output$M3sig_5_taxonlist <- renderText({
          
          sig.tab5 <- kable(tab.five, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab5
        })
        
        
        incProgress(1/10, message = "Displaying Results in progress")
        output$taxa_surv_downloadTable = renderUI({
          tagList(
            p(" ", style = "margin-top: 20px;"),
            box(title = strong("Download Output Table"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"),
                #h5("Summary Statistics"),
                #downloadButton("tdownloadTabl1", "Download", width = '50%', style = "background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("tdownloadTabl2", "Download", width = '50%', style = "background-color: red3"),
                h5("Dendrogram"),
                downloadButton("gdownload", "Download", width = '50%', style = "background-color: red3")
            )
          )
        })
        
        output$tdownloadTabl2 <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Output.zip")
          },
          content = function(DA.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.cox.out$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.cox.out$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.cox.out$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.cox.out$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.cox.out$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.cox.out$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=DA.file, files=dataFiles)
          }
        )
        
        output$gdownload <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Graphical.Output", ".html", sep="")
          },
          content = function(file) {
            htmlwidgets::saveWidget(as_widget(taxa.sig.dend(taxa.cox.out, chooseData$NAadded$tax.tab, "twopi", include)), file)
          }
        )
        
        ref_string = REFERENCE_CHECK(data_transform = input$surv.dataType_taxa, method_name = isolate(input$surv3.taxa.method), FDR = "Yes")
        if (is.null(ref_string)) {
          shinyjs::hide("referencesM3.taxa")
        } else {
          shinyjs::show("referencesM3.taxa")
          output$referencesM3.taxa = renderUI({
            tagList(
              box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
      }
    )
    
    
    shinyjs::enable("runbtn_CoxT")
  })
  
  
  ###################################################
  ####### Taxa Random Survival Forest Analysis ######
  ###################################################
  observeEvent(input$runbtn_model4, {
    
    taxa.sam.dat <- chooseData$sam.dat[match(rownames(chooseData$taxa.out$clr$phylum),rownames(chooseData$sam.dat)),]
    surv.dat <- surv.events.dat.func(taxa.sam.dat, 
                                     surv.con = chooseData$surv.Time.select, 
                                     surv.second.bin = chooseData$censor.select)
    taxa.dat.4 <- chooseData$taxa.out[[taxa.types$dataType_model4]]
    
    if (input$subgroup.4 == "Yes"){
      
      ind.sub <- taxa.sam.dat[[input$subgroup.sel.4]] == input$pick_subgroup.4
      taxa.sam.dat <- taxa.sam.dat[-ind.sub,]  
      surv.dat <- surv.events.dat.func(taxa.sam.dat, 
                                       surv.con = chooseData$surv.Time.select, 
                                       surv.second.bin = chooseData$censor.select)
      
      taxa.dat.4 <- lapply(taxa.dat.4, function(x){
        return(x[-ind.sub,])})
    }
    
    validate(
      if (sum(surv.dat$status, na.rm = T) <= 0) {
        showNotification("Error: No Events in Data",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtn_model4")
    shinyjs::disable("surv4.Time.select")
    shinyjs::disable("surv4.dataType_taxa")
    shinyjs::disable("surv4.censor.select")
    shinyjs::disable("surv4.method.select")
    shinyjs::disable("surv4.fold")
    
    withProgress(
      
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        incProgress(1/10, message = "Calculating")
        
        rank <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
        
        out <- list()
        
        folds <- as.numeric(isolate(input$surv4.fold.select))
        trees <- as.numeric(isolate(input$surv4.tree.select))
        
        #taxa.dat <- chooseData$taxa.out[[taxa.types$dataType_model4]]
        
        if(input$include_speciesM4 == "Phylum - Genus (Default)"){
          ranks = 5
        } else {
          ranks = 6
        }
        
        for(i in 1:ranks) {
          taxa.dat.4[[i]] <- taxa.dat.4[[i]][rownames(surv.dat),]
        }
        
        
        if (input$surv4.method.select == "Random survival forests") {
          incProgress(6/10, message = "Random survival forests")
          set.seed(521)
          rf.fit.out <- try(surv.random.forest(taxa.dat.4, chooseData$taxa.names.out, surv.dat, n.tree = trees, ranks.upto = ranks), silent = TRUE)
          
          importance <- list()
          duplicates <- list()
          err.dat <- list()
          
          for (k in 1:ranks) {
            imp <- as.matrix(rf.fit.out$rf.fit[[k]]$importance)
            colnames(imp) <- "Importance"
            importance[[k]] <- imp
            
              
            
            if (sum(is.na(chooseData$taxa.names.out$duplicates[[k]])) == 0) {
              duplicate.taxa <- sapply(strsplit(unlist(chooseData$taxa.names.out$duplicates[[k]]), " :"),  "[", 1)
              duplicates[[k]] <- rep(TRUE, length(duplicate.taxa)) #%in% taxon.inplot
            } else {
              duplicates[[k]] <- FALSE
            }
            err.dat[[k]] = try(as.data.frame(cbind(ntree = c(1:trees)[!is.na(rf.fit.out$rf.fit[[k]]$err.rate)], err.rate = rf.fit.out$rf.fit[[k]]$err.rate[!is.na(rf.fit.out$rf.fit[[k]]$err.rate)])), silent = TRUE)
          }
          incProgress(2/10, message = "Plotting")
          output$surv4_display_results = renderUI({
            tagList(
              do.call(tabsetPanel, lapply(1:ranks, function(i) {
                tabPanel(title = rank[i], align = "center",
                         tabsetPanel(tabPanel(title = "Variable Importance", align = "center",
                                              plotOutput(paste0("importance", i), height = 600, width = 500),
                                              plotOutput(paste0("duplicates_m4", i), height = 11.5*sum(duplicates[[i]])+10, width = 500)),
                                     tabPanel(title = "Error Rate", align = "center",
                                              plotOutput(paste0("err.rate", i), height = 600, width = 500)))
                )
              }))
            )
          })
          
          
          incProgress(1/10, message = "Plotting")
          lapply(1:ranks, function(j) {
            if (class(rf.fit.out)[1] == "try-error") {
              output[[paste0("importance", j)]] <- renderPlot({
                text(x = 0.5, y = 0.5, "Error : Please try different number of folds or follow-up period.", 
                     cex = 1.2, col = "black")
              })
            } else {
              output[[paste0("importance", j)]] <- renderPlot({
                plot.importance(rf.fit.out$rf.fit[[j]], input$surv4.num.display)
              })
              output[[paste0("err.rate", j)]] <- renderPlot({
                par(mar = c(7,9,2,2))
                plot(smooth.spline(err.dat[[j]]), type = "l", ylab = "Error Rate", xlab = "Number of Trees", lwd = 1.5)
              })
              
              lapply(1:ranks, function(j) {
                output[[paste0("duplicates_m4", j)]] <- renderPlot({
                  surv.en.duplicate.list(duplicates[[j]], chooseData$taxa.names.out$duplicates[[j]])
                })
              })
              
              output$downloadTable_m4 = renderUI({
                tagList(
                  p(" ", style = "margin-top: 20px;"),
                  box(title = strong("Download Output Table"), width = NULL, status = "primary", solidHeader = TRUE,
                      p("You can download the data analysis outputs.",
                        style = "font-size:11pt"),
                      h5("Data Analysis Outputs"),
                      downloadButton("m4tdownloadTabl", "Download", width = '50%', style = "background-color: red3")
                  )
                )
              })
              
              output$m4tdownloadTabl <- downloadHandler(
                filename = function() {
                  paste("Model4.Analysis.Output.zip")
                },
                content = function(DA.file) {
                  temp <- setwd(tempdir())
                  on.exit(setwd(temp))
                  if (ranks == 5) {
                    dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
                    write.table(as.data.frame(importance[[1]]), file = "Phylum.txt", sep = "\t")
                    write.table(as.data.frame(importance[[2]]), file = "Class.txt", sep = "\t")
                    write.table(as.data.frame(importance[[3]]), file = "Order.txt", sep = "\t")
                    write.table(as.data.frame(importance[[4]]), file = "Family.txt", sep = "\t")
                    write.table(as.data.frame(importance[[5]]), file = "Genus.txt", sep = "\t")
                  }
                  else if(ranks == 6) {
                    dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
                    write.table(as.data.frame(importance[[1]]), file = "Phylum.txt", sep = "\t")
                    write.table(as.data.frame(importance[[2]]), file = "Class.txt", sep = "\t")
                    write.table(as.data.frame(importance[[3]]), file = "Order.txt", sep = "\t")
                    write.table(as.data.frame(importance[[4]]), file = "Family.txt", sep = "\t")
                    write.table(as.data.frame(importance[[5]]), file = "Genus.txt", sep = "\t")
                    write.table(as.data.frame(importance[[6]]), file = "Species.txt", sep = "\t")
                  }
                  zip(zipfile=DA.file, files=dataFiles)
                }
              )
            }
          })
        } 
        
        ref_string = REFERENCE_CHECK(data_transform = input$surv4.dataType_taxa, method_name = isolate(input$surv4.method.select), FDR = "No")
        if (is.null(ref_string)) {
          shinyjs::hide("referencesM4")
        } else {
          shinyjs::show("referencesM4")
          output$referencesM4 = renderUI({
            tagList(
              box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
      }
    )
    
    
    shinyjs::enable("surv4.censor.select")
    shinyjs::enable("surv4.method.select")
    shinyjs::enable("surv4.fold.select")
    shinyjs::enable("surv4.Time.select")
    shinyjs::enable("surv4.dataType_taxa")
    shinyjs::enable("runbtn_model4")
  })
  
  
  ##########################################
  ##                  QC                  ##
  ##########################################
  observeEvent(input$run, {
    
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        
        incProgress(1/10, message = "Data Trimming in progress")
        
        if (nchar(input$part.rem.str) == 0) {
          rem.tax.complete <- rem.tax.d
          rem.tax.partial <- rem.tax.str.d
        } else {
          rem.tax.complete <- unique(c(unlist(strsplit(input$rem.str, split = ",")), rem.tax.d))
          rem.tax.partial <- unique(c(unlist(strsplit(input$part.rem.str, split = ",")), rem.tax.str.d))
        }
        
        tax.tab <- tax_table(infile$biom)
        
        if (input$kingdom != "all") {
          ind <- is.element(tax.tab[,1], input$kingdom)
          validate(
            if (sum(ind) == 0) {
              showNotification(h4(paste("Error: Please select valid Kingdom. Available kingdoms are:", 
                                        paste(c(na.omit(unique(tax.tab[,1])) ,"and all"), collapse = ", "))),
                               type = "error")
            } else {
              NULL
            }
          )
        }
        
        shinyjs::disable("run")
        shinyjs::disable("slider1")
        shinyjs::disable("slider2")
        shinyjs::disable("kingdom")
        shinyjs::disable("skip")
        shinyjs::disable("binwidth")
        shinyjs::disable("binwidth2")
        
        rcol$selected = "rgba(255, 0, 0, 0.6)"
        
        chooseData$surv.Time.select = input$surv.Time.select
        #chooseData$subgroup = input$subgroup
        chooseData$censor.select = input$censor.select
        chooseData$follow.time = input$follow.time
        
        infile$qc_biom = biom.cleanS(infile$biom,
                                     input$kingdom,
                                     lib.size.cut.off = input$slider1,
                                     mean.prop.cut.off = input$slider2/100,
                                     rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
                                     surv.time = chooseData$surv.Time.select,
                                     censor = chooseData$censor.select,
                                     follow = chooseData$follow.time)
        infile$qc_biomNA = biom.cleanSNA(infile$biom,
                                         input$kingdom,
                                         lib.size.cut.off = input$slider1,
                                         mean.prop.cut.off = input$slider2/100,
                                         rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
                                         surv.time = chooseData$surv.Time.select,
                                         censor = chooseData$censor.select,
                                         follow = chooseData$follow.time)
        
        # if (input$subgroup == "No") {
        #   infile$qc_biom = biom.cleanS(infile$biom,
        #                                input$kingdom,
        #                                lib.size.cut.off = input$slider1,
        #                                mean.prop.cut.off = input$slider2/100,
        #                                rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
        #                                surv.time = chooseData$surv.Time.select,
        #                                censor = chooseData$censor.select,
        #                                follow = chooseData$follow.time)
        #   infile$qc_biomNA = biom.cleanSNA(infile$biom,
        #                                    input$kingdom,
        #                                    lib.size.cut.off = input$slider1,
        #                                    mean.prop.cut.off = input$slider2/100,
        #                                    rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
        #                                    surv.time = chooseData$surv.Time.select,
        #                                    censor = chooseData$censor.select,
        #                                    follow = chooseData$follow.time)
        # }
        # else if (input$subgroup == "Yes") {
        #   chooseData$subgroup.sel = input$subgroup.sel
        #   chooseData$pick_subgroup = input$pick_subgroup
        # 
        #   infile$qc_biom = biom.cleanS(infile$biom,
        #                                input$kingdom,
        #                                lib.size.cut.off = input$slider1,
        #                                mean.prop.cut.off = input$slider2/100,
        #                                rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
        #                                surv.time = chooseData$surv.Time.select,
        #                                censor = chooseData$censor.select,
        #                                follow = chooseData$follow.time,
        #                                subgroup.var = chooseData$subgroup.sel,
        #                                subgroup = chooseData$pick_subgroup)
        #   infile$qc_biomNA = biom.cleanSNA(infile$biom,
        #                                input$kingdom,
        #                                lib.size.cut.off = input$slider1,
        #                                mean.prop.cut.off = input$slider2/100,
        #                                rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
        #                                surv.time = chooseData$surv.Time.select,
        #                                censor = chooseData$censor.select,
        #                                follow = chooseData$follow.time,
        #                                subgroup.var = chooseData$subgroup.sel,
        #                                subgroup = chooseData$pick_subgroup)
        # }
        
        
        incProgress(3/10, message = "Rarefying in progress")
        lib_size.sum = lib.size.func(infile$qc_biom)$lib.size.sum
        infile$rare_biom = rarefy.func(infile$qc_biom, 
                                       cut.off = lib_size.sum["Minimum"],
                                       multi.rarefy = 1)
        infile$rare_biomNA = rarefy.func(infile$qc_biomNA, 
                                         cut.off = lib_size.sum["Minimum"],
                                         multi.rarefy = 1)
        
        incProgress(2/10, message = "Saving File in progress")
        
        chooseData$sam.dat = sample_data(infile$qc_biom)
        chooseData$mon.sin.rev.bin.con = is.mon.sin.rev.bin.con(chooseData$sam.dat)
        chooseData$prim_vars = surv.pri.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con)$is.bin
        chooseData$prim_varsSurv = surv.pri.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con)
        chooseData$tax.tab = tax_table(infile$rare_biom)
        chooseData$tax.tabNA = tax_table(infile$rare_biomNA)
        chooseData$prim_varsCox = surv.pri.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con)
        
        chooseData$nt.selected.bin <- chooseData$prim_varsSurv[[1]][ !chooseData$prim_varsSurv[[1]] %in% c(input$censor.select)]
        chooseData$nt.selected.con <- chooseData$prim_varsSurv[[2]][ !chooseData$prim_varsSurv[[2]] %in% c(input$surv.Time.select)]
        
        
        #library.size <- library.size[names(library.size) %in% rownames(rare.sam.dat)]
        
        output$moreControls <- renderUI({
          tagList(
            box(title = strong("Download Data"), width = NULL, status = "primary", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:13pt"),
                h5("Data after Quality Control"),
                downloadButton("downloadData2", "Download", width = '50%', style = "background-color: red3"),br(),
                h5("Data after Quality Control and Rarefaction"),
                downloadButton("downloadData3", "Download", width = '50%', style = "background-color: red3"),br(),
                p("For your reference, you can download the data files above for the phyloseq object (biom.after.qc) after QC and
                      (rare.biom.after.qc) after QC and rarefaction.",
                  style = "font-size:11pt")
            ), 
            uiOutput("ref_micloud_qc")
          )
        })
        
        output$text <- renderText({"You are all set! You can proceed to data analysis!"})
        output$ref_micloud_qc <- renderUI({
          tagList(
            box(title = strong("Reference"), width = NULL, status = "primary", solidHeader = TRUE, 
                p("Gu W, Moon J, Chisina C, Kang B, Park T, Koh H. MiCloud: A unified web platform for comprehensive microbiome data analysis. PLoS One 2022:17(8): e0272354.")
            )
          )
        })
        
        qc_biom = infile$qc_biom
        output$downloadData2 <- downloadHandler(
          filename = function() {
            paste("biom.after.qc.Rdata")
          },
          content = function(file1) {
            save(qc_biom, file = file1)
          })
        
        rare_biom = infile$rare_biom
        output$downloadData3 <- downloadHandler(
          filename = function() {
            paste("rare.biom.after.qc.Rdata")
          },
          content = function(file1) {
            save(rare_biom, file = file1)
          })
        
        incProgress(1/10, message = "Done")
        shinyjs::enable("run")
        shinyjs::enable("slider1")
        shinyjs::enable("slider2")
        shinyjs::enable("kingdom")
        shinyjs::enable("skip")
        shinyjs::enable("binwidth")
        shinyjs::enable("binwidth2")
      })
    
    
    # shinyjs::disable("covariatesCoxA")
    # shinyjs::disable("surv3.alpha.method")
    shinyjs::disable("runbtn_CoxA")
    # shinyjs::disable("chooseAdjustment3A")
    # shinyjs::disable("covariatesOptionsCoxA")
    
    
    # shinyjs::disable("covariatesCoxB")
    # shinyjs::disable("surv3.beta.method")
    shinyjs::disable("runbtn_CoxB")
    # shinyjs::disable("covariatesOptionsCoxB")
    
    # shinyjs::disable("primvar")
    # shinyjs::disable("chooseMethod")
    shinyjs::disable("runbtn_bin")
    # shinyjs::disable("covariates")
    # shinyjs::disable("chooseAdjustment")
    
    # shinyjs::disable("beta_chooseMethod_bin")
    shinyjs::disable("beta_runbtn_cross_bin")
    # shinyjs::disable("beta.primvar_cross")
    # shinyjs::disable("beta_covariates_bin")
    
    # shinyjs::disable("chooseMethod_taxa")
    shinyjs::disable("taxa_runbtn_bin")
    # shinyjs::disable("surv3.taxa.method")
    shinyjs::disable("runbtn_CoxT")
    # shinyjs::disable("surv4.method.select")
    shinyjs::disable("runbtn_model4")
  })
  
  ##########################################
  ## Data Calculation (Alpha, Beta, Taxa) ##
  ##########################################
  observeEvent(input$divCalcRun, {
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        shinyjs::disable("divCalcRun")
        
        rare.otu.tab <- otu_table(infile$rare_biom)
        
        incProgress(3/10, message = "Calculating Diversity")
        
        chooseData$alpha.div.rare = alpha.v1.func(infile$rare_biom)
        chooseData$alpha.div.qc = alpha.v1.func(infile$qc_biom)
        chooseData$alpha.div = chooseData$alpha.div.rare
        #chooseData$alpha.div <- chooseData$alpha.div[rownames(chooseData$alpha.div) %in% rownames(chooseData$sam.dat),]
        
        incProgress(3/10, message = "Calculating Distance")
        
        ds.Ks$res = Ds.Ks.func(infile$rare_biom, infile$qc_biom)
        
        incProgress(1/10, message = "Transforming Data")
        
        rare.otu.tab <- otu_table(infile$rare_biom)
        rare.tax.tab <- tax_table(infile$rare_biom)
        no.rare.otu.tab <- otu_table(infile$qc_biom)
        no.rare.tax.tab <- tax_table(infile$qc_biom)
        
        rare.otu.tabNA <- otu_table(infile$rare_biomNA)
        rare.tax.tabNA <- tax_table(infile$rare_biomNA)
        no.rare.otu.tabNA <- otu_table(infile$qc_biomNA)
        no.rare.tax.tabNA <- tax_table(infile$qc_biomNA)
        
        chooseData$taxa.out = tax.trans(no.rare.otu.tab, no.rare.tax.tab, rare.otu.tab, rare.tax.tab)
        chooseData$taxa.outNA = tax.trans.na(no.rare.otu.tabNA, no.rare.tax.tabNA, rare.otu.tabNA, rare.tax.tabNA)
        chooseData$taxa.names.out = taxa.names.rank(chooseData$taxa.out[[1]])
        chooseData$tax.tab = rare.tax.tab
        chooseData$tax.tabNA = rare.tax.tabNA
        chooseData$NAadded <- add_NA(chooseData$taxa.outNA, chooseData$tax.tabNA)
        
        output$divCalcDownload <- renderUI({
          tagList(
            box(title = strong("Download Data"), width = NULL, status = "primary", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:15pt"),
                p("You can download alpha- and beta-diversity data.",
                  style = "font-size:11pt"), 
                h5("Alpha Diversity", HTML('&emsp;'), HTML('&nbsp;'), HTML('&nbsp;'), "Beta Diversity"),
                downloadButton("alphaDiv", "Download", width = '50%', style = "background-color: red3"), HTML('&emsp;'),
                downloadButton("betaDiv", "Download", width = '50%', style = "background-color: red3"),br(),br(),
                p("You can download taxonomic abundance data.",
                  style = "font-size:11pt"),
                h5("Count", HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), "Count (Rarefied)"),
                downloadButton("taxadataCount", "Download", width = '50%', style = "background-color: red3"), HTML('&emsp;'),
                downloadButton("taxadataRareCount", "Download", width = '50%', style = "background-color: red3"), br(), 
                h5("Proportion", HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&nbsp;'), "CLR"),
                downloadButton("taxadataProp", "Download", width = '50%', style = "background-color: red3"), HTML('&emsp;'),
                downloadButton("taxadataCLR", "Download", width = '50%', style = "background-color: red3"), br(),
                h5("Arcsine-root"),
                downloadButton("taxadataArc", "Download", width = '50%', style = "background-color: red3"), br(),br(),
            )
          )
        })
        
        alpha.div = chooseData$alpha.div
        
        output$alphaDiv <- downloadHandler(
          filename = function() {
            paste("Alpha.Diversity.txt")
          },
          content = function(alpha.file) {
            write.table(chooseData$alpha.div, file = alpha.file, row.names = TRUE, col.names = TRUE, sep = "\t")
          })
        
        output$betaDiv <- downloadHandler(
          filename = function() {
            paste("Beta.Diversity.zip")
          },
          content <- function(fname) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Jaccard.txt", "Bray.Curtis.txt", "U.UniFrac.txt" ,"G.UniFrac.txt", "W.UniFrac.txt")
            
            write.table(as.data.frame(ds.Ks$res$Ds$Jaccard), file = "Jaccard.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$Bray.Curtis), file = "Bray.Curtis.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$U.UniFrac), file = "U.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$G.UniFrac), file = "G.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$W.UniFrac), file = "W.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=fname, files=dataFiles)
          }
        )
        
        
        output$datTransDownload <- renderUI({
          tagList(
            box(title = strong("Download Data"), width = NULL, status = "primary", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:15pt"),
                p("You can download taxonomic abundance data.",
                  style = "font-size:11pt"), 
                h5("Count"),
                downloadButton("taxadataCount", "Download", width = '50%', style = "background-color: red3"), br(), 
                h5("Count (Rarefied)"),
                downloadButton("taxadataRareCount", "Download", width = '50%', style = "background-color: red3"), br(), 
                h5("Proportion"),
                downloadButton("taxadataProp", "Download", width = '50%', style = "background-color: red3"), br(), 
                h5("CLR"),
                downloadButton("taxadataCLR", "Download", width = '50%', style = "background-color: red3"), br(),
                h5("Arcsine-root"),
                downloadButton("taxadataArc", "Download", width = '50%', style = "background-color: red3"), br(),br(),
            )
          )
        })
        
        count_biom = chooseData$taxa.out$count
        rare_biom = chooseData$taxa.out$rare.count
        prop_biom = chooseData$taxa.out$prop
        clr_biom = chooseData$taxa.out$clr
        arc_biom = chooseData$taxa.out$arcsin
        
        output$taxadataCount <- downloadHandler(
          
          filename = function() {
            paste("Count.Data.zip")
          },
          content = function(count.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(count_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=count.file, files=dataFiles)
          }
        )
        
        incProgress(3/10, message = "Saving")
        
        output$taxadataRareCount <- downloadHandler(
          
          filename = function() {
            paste("Rarefied.Count.Data.zip")
          },
          content = function(rare.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(rare_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=rare.file, files=dataFiles)
          }
        )
        output$taxadataProp <- downloadHandler(
          filename = function() {
            paste("Proportion.Data.zip")
          },
          content = function(prop.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(prop_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=prop.file, files=dataFiles)
          }
        )
        output$taxadataCLR <- downloadHandler(
          filename = function() {
            paste("CLR.Transformed.Data.zip")
          },
          content = function(clr.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(clr_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=clr.file, files=dataFiles)
          }
        )
        
        output$taxadataArc <- downloadHandler(
          filename = function() {
            paste("Arcsin.Transformed.Data.zip")
          },
          content = function(clr.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(arc_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=clr.file, files=dataFiles)
          }
        )
        
        incProgress(1/10, message = "Done")
        shinyjs::enable("divCalcRun")
      })
    
    shinyjs::enable("covariatesCoxA")
    shinyjs::enable("surv3.alpha.method")
    shinyjs::enable("runbtn_CoxA")
    shinyjs::enable("covariatesCoxB")
    shinyjs::enable("surv3.beta.method")
    shinyjs::enable("runbtn_CoxB")
    shinyjs::enable("primvar")
    shinyjs::enable("chooseMethod")
    shinyjs::enable("runbtn_bin")
    shinyjs::enable("beta_chooseMethod_bin")
    shinyjs::enable("beta_runbtn_cross_bin")
    shinyjs::enable("chooseMethod_taxa")
    shinyjs::enable("taxa_runbtn_bin")
    shinyjs::enable("surv3.taxa.method")
    shinyjs::enable("runbtn_CoxT")
    shinyjs::enable("surv4.method.select")
    shinyjs::enable("runbtn_model4")
  })
  
  ##########################################
  ## Alpha Cross-Sectional Data Analysis ##
  ##########################################
  observeEvent(input$runbtn_bin,{
    validate(
      if (input$covariates == "Covariate(s)" & is.null(input$covariatesOptions)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtn_bin")
    shinyjs::disable("chooseAdjustment")
    shinyjs::disable("primvar")
    shinyjs::disable("chooseMethod")
    shinyjs::disable("covariates")
    shinyjs::disable("covariatesOptions")
    shinyjs::disable("alphaCat1")
    shinyjs::disable("alphaCat2")
    shinyjs::disable("chooseData")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        incProgress(1/10, message = "Renaming")
        chooseData$alpha.div = chooseData$alpha.div.rare
        
        alpha.categors = c(alpha.categos$cat1, alpha.categos$cat2)
        alpha.bin_categos = c(input$alphaCat1, input$alphaCat2)
        
        rename.cats_ref = alpha.bin_categos[which(alpha.categors == alpha.categos$cat1)]
        rename.cats_com = alpha.bin_categos[which(alpha.categors != alpha.categos$cat1)]
        
        alpha.bin.cat.ref.ori.out <- alpha.bin.cat.ref.ori.func(chooseData$sam.dat, input$primvar)
        sam_dat <- alpha.bin.cat.recode.func(chooseData$sam.dat, input$primvar, alpha.bin.cat.ref.ori.out,
                                             rename.cats_ref, rename.cats_com)
        
        incProgress(3/10, message = "Calculating")
        alpha.div.bin.out <- alpha.bin.cat.ref.func(input$primvar, rename.cats_ref,
                                                    rename.cats_com, sam_dat,
                                                    chooseData$alpha.div)
        
        alpha.results$bin.var <- alpha.div.bin.out$bin.var
        alpha.results$alpha_div <- alpha.div.bin.out$alpha.div
        alpha.results$alpha.bin.sum.out <- alpha.bin.sum.func(bin.var = alpha.div.bin.out$bin.var,
                                                              alpha.div = alpha.div.bin.out$alpha.div)
        
        if (input$chooseMethod == "Welch t-test" | input$chooseMethod == "Wilcoxon rank-sum test (default)") {
          if (input$chooseMethod == "Welch t-test") {
            incProgress(3/10, message = "T test")
            t.test.out <- alpha.bin.t.test(alpha.results$bin.var, alpha.results$alpha_div)
            alpha.data.results$table.p.out = t.test.out
            alpha.data.results$data.q.out = q.func(t.test.out, method = "BH")
          }
          else if (input$chooseMethod == "Wilcoxon rank-sum test (default)") {
            incProgress(3/10, message = "Wilcoxon test")
            wilcox.test.out <- alpha.bin.wilcox.test(alpha.results$bin.var, alpha.results$alpha_div)
            alpha.data.results$table.p.out = wilcox.test.out
            alpha.data.results$data.q.out = q.func(wilcox.test.out, method = "BH")
          }
          
          multi.test$boolval = FALSE
          alpha.data.results$table.out = alpha.data.results$table.p.out
          
          # if (input$chooseAdjustment == "Yes") {
          #   multi.test$boolval = TRUE
          #   alpha.data.results$table.out = alpha.data.results$data.q.out
          # }
          # else if (input$chooseAdjustment == "No (Default)") {
          #   multi.test$boolval = FALSE
          #   alpha.data.results$table.out = alpha.data.results$table.p.out
          # }
          
          incProgress(3/10, message = "Graph")
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong("Graphs for Alpha Diversity"), 
                  align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                  plotOutput("box_plots", height = 850, width = 500)
              )
            )
          })
          output$box_plots = renderPlot({
            alpha.bin.hist(alpha.results$bin.var, alpha.results$alpha_div, alpha.data.results$data.q.out, multi.test$boolval)
          })
        }
        else if (input$chooseMethod == "Linear regression") {
          alpha.bin.cov.out <- alpha.bin.cov.cat.ref.func(input$primvar, rename.cats_ref,
                                                          rename.cats_com, input$covariatesOptions,
                                                          sam_dat, chooseData$alpha.div)
          alpha.reg.results$bin.var <- alpha.bin.cov.out$bin.var
          alpha.reg.results$cov.var <- alpha.bin.cov.out$cov.var
          alpha.reg.results$alpha.div <- alpha.bin.cov.out$alpha.div
          
          if (input$chooseMethod == "Linear regression") {
            incProgress(3/10, message = "Linear regression with Covariate(s)")
            
            alpha.lm.bin.cov.out <- alpha.lm.bin.cov.func(bin.var = alpha.reg.results$bin.var,
                                                          cov.var = alpha.reg.results$cov.var,
                                                          alpha.div = alpha.reg.results$alpha.div,
                                                          scale = TRUE)
            alpha.data.results$table.p.out = alpha.lm.bin.cov.out
            alpha.data.results$data.q.out = q.func(alpha.lm.bin.cov.out, method = "BH")
            alpha.data.results$table.out = alpha.data.results$data.q.out
          }
          
          multi.test$boolval = TRUE
          alpha.data.results$table.out = alpha.data.results$data.q.out
          
          # if (input$chooseAdjustment == "Yes") {
          #   multi.test$boolval = TRUE
          #   alpha.data.results$table.out = alpha.data.results$data.q.out
          # }
          # else if (input$chooseAdjustment == "No (Default)") {
          #   multi.test$boolval = FALSE
          #   alpha.data.results$table.out = alpha.data.results$table.p.out
          # }
          
          incProgress(3/10, message = "Graph")
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong("Graphs for Alpha Diversity"), 
                  align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                  plotOutput("forest_plots", height = 400, width = 500)
              )
            )
          })
          
          output$forest_plots = renderPlot({
            alpha.forest.plot(alpha.data.results$data.q.out, multi.test$boolval)
          })
          
        }
        
        output$alpha_downloadTable = renderUI({
          tagList(
            p(" ", style = "margin-top: 20px;"),
            box(title = strong("Download Output Table"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Summary Statistics"),
                downloadButton("downloadTabl1", "Download", width = '50%', style = "background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("downloadTabl2", "Download", width = '50%', style = "background-color: red3")
            )
          )
        })
        output$downloadTabl1 <- downloadHandler(
          filename = function() {
            paste("Alpha.Sum.Table.txt")
          },
          content = function(file) {
            write.table(alpha.results$alpha.bin.sum.out, file, sep="\t")
          }
        )
        output$downloadTabl2 <- downloadHandler(
          filename = function() {
            paste("Alpha.DA.Output.txt")
          },
          content = function(file) {
            write.table(alpha.data.results$table.out, file)
          }
        )
        
        ref_string = REFERENCE_CHECK(method_name = isolate(input$chooseMethod), FDR = "No (Default)")
        if (is.null(ref_string)) {
          shinyjs::hide("referencesM2.alpha")
        } else {
          shinyjs::show("referencesM2.alpha")
          output$referencesM2.alpha = renderUI({
            tagList(
              box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        
        shinyjs::enable("runbtn_bin")
        shinyjs::enable("primvar")
        shinyjs::enable("chooseAdjustment")
        shinyjs::enable("chooseMethod")
        shinyjs::enable("covariates")
        shinyjs::enable("covariatesOptions")
        shinyjs::enable("alphaCat1")
        shinyjs::enable("alphaCat2")
        shinyjs::enable("chooseData")
      })
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  ##########################################
  ## Beta Cross-Sectional Data Analysis ###
  ##########################################
  observeEvent(input$beta_runbtn_cross_bin,{
    validate(
      if (input$beta_covariates_bin == "Covariate(s)" & is.null(input$beta_covariatesOptions_bin)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("beta_runbtn_cross_bin")
    shinyjs::disable("beta.primvar_cross")
    shinyjs::disable("betaCat1")
    shinyjs::disable("betaCat2")
    shinyjs::disable("beta_covariates_bin")
    shinyjs::disable("beta_covariatesOptions_bin")
    shinyjs::disable("beta_chooseMethod_bin")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        incProgress(3/10, message = "Rename")
        
        beta.categors = c(beta.categos$cat1, beta.categos$cat2)
        beta.bin_categos = c(input$betaCat1, input$betaCat2)
        
        rename.catsbin_ref = beta.bin_categos[which(beta.categors == beta.categos$cat1)]
        rename.catsbin_com = beta.bin_categos[which(beta.categors != beta.categos$cat1)]
        
        beta.bin.cat.ref.ori.out <- beta.bin.cat.ref.ori.func(chooseData$sam.dat, input$beta.primvar_cross)
        beta.sam_dat.bin <- beta.bin.cat.recode.func(chooseData$sam.dat, input$beta.primvar_cross,
                                                     beta.bin.cat.ref.ori.out,
                                                     rename.catsbin_ref, rename.catsbin_com)
        
        beta.sam_dat.bin <- sam_no_miss_cov(beta.sam_dat.bin, input$beta_covariatesOptions_bin)
        
        if (input$beta_covariates_bin == "None") {
          
          if (input$beta_chooseMethod_bin == "MiRKAT") {
            
            incProgress(3/10, message = "MiRKAT without Covariate(s)")
            beta.bin.out <- beta.bin.cat.ref.func(input$beta.primvar_cross,
                                                  rename.catsbin_ref, rename.catsbin_com,
                                                  beta.sam_dat.bin, Ds.Ks = ds.Ks$res)
            beta.data.results$data.q.out <- beta.bin.out
          }
        } else if (input$beta_covariates_bin == "Covariate(s)") {
          if (is.null(input$beta_covariatesOptions_bin)) {
            if (input$beta_chooseMethod_bin == "MiRKAT") {
              incProgress(3/10, message = "MiRKAT without Covariate(s)")
              beta.bin.out <- beta.bin.cat.ref.func(input$beta.primvar_cross,
                                                    rename.catsbin_ref, rename.catsbin_com,
                                                    beta.sam_dat.bin, Ds.Ks = ds.Ks$res)
              beta.data.results$data.q.out <- beta.bin.out
            }
          } else {
            
            beta_re_cov <- beta_no_miss_cov(beta.sam_dat.bin, ds.Ks$res,input$beta_covariatesOptions_bin)
            
            beta.bin.cov.out <- beta.bin.cov.cat.ref.func(input$beta.primvar_cross,
                                                          rename.catsbin_ref, rename.catsbin_com,
                                                          input$beta_covariatesOptions_bin, beta.sam_dat.bin,
                                                          Ds.Ks = beta_re_cov)
            
            
            if (input$beta_chooseMethod_bin == "MiRKAT") {
              incProgress(3/10, message = "MiRKAT with Covariate(s)")
              beta.data.results$data.q.out <- beta.bin.cov.out
            }
          }
        }
        
        incProgress(3/10, message = "Plotting Graph")
        
        if (input$beta_covariates_bin == "None") {
          beta.down.results$CS = mirkat.bin(beta.data.results$data.q.out)
          beta.plot.info <- mirkat.bin.plot1(beta.down.results$CS, beta.data.results$data.q.out)
          output$m2_beta_2d = renderPlot({
            isolate(try(mirkat.bin.plot2(beta.down.results$CS, beta.data.results$data.q.out, beta.plot.info$mod, beta.plot.info$sub.tit), silent = TRUE))
          })
        } else if (input$beta_covariates_bin == "Covariate(s)") {
          beta.down.results$CS = mirkat.bin.cov(beta.data.results$data.q.out)
          beta.plot.info <- mirkat.bin.plot1(beta.down.results$CS, beta.data.results$data.q.out)
          output$m2_beta_2d = renderPlot({
            isolate(try(mirkat.bin.plot2(beta.down.results$CS, beta.data.results$data.q.out, beta.plot.info$mod, beta.plot.info$sub.tit), silent = TRUE))
          })
        }
        
        beta.data.results$beta.plot.info <- beta.plot.info
        
        output$beta_display_results_cross = renderUI({
          tagList(
            box(title = strong("Graphs for Beta Diversity"), 
                align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("m2_beta_2d", height = 850, width = 500)
            )
          )
        })
        
        beta_ind <- c("Jaccard", "Bray-Curtis", "U.UniFrac", "G.UniFrac", "W.UniFrac")
        
        output$beta_downloadTable = renderUI({
          tagList(
            p(" ", style = "margin-top: 20px;"),
            box(title = strong("Download Output Table"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Data Analysis Outputs"),
                downloadButton("beta_downloadTabl1", "Download", width = '50%', style = "background-color: red3"),
                h5("Data Analysis Plot Image"),
                downloadButton("beta_downloadTabl2", "Download", width = '50%', style = "background-color: red3")
            )
          )
        })
        incProgress(3/10, message = "SAVE")
        
        output$beta_downloadTabl1 <- downloadHandler(
          filename = function() {
            paste("Beta.DA.Output.txt")
          },
          content = function(file) {
            out_temp = as.data.frame(unlist(beta.down.results$CS))
            rownames(out_temp) = c("Jaccard","Bray.Curtis","U.UniFrac","G.UniFrac","W.UniFrac", "omnibus_p")
            colnames(out_temp) = "p-value"
            beta.down.results$CS = out_temp
            write.table(beta.down.results$CS, file, sep="\t")
          }
        )
        
        output$beta_downloadTabl2 <- downloadHandler(
          filename = function() {
            paste("Beta.DA.Output.png")
          },
          content = function(file) {
            png(file=file,
                width=1050, height=950)
            print( mirkat.bin.plot2(beta.down.results$CS, beta.data.results$data.q.out, beta.plot.info$mod, beta.plot.info$sub.tit) )
            dev.off()
          }
        )
        
        ref_string = REFERENCE_CHECK(method_name = isolate(input$beta_chooseMethod_bin))
        if (is.null(ref_string)) {
          shinyjs::hide("referencesM2.beta")
        } else {
          shinyjs::show("referencesM2.beta")
          output$referencesM2.beta = renderUI({
            tagList(
              box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
      }
    )
    
    shinyjs::enable("beta_runbtn_cross_bin")
    shinyjs::enable("beta.primvar_cross")
    shinyjs::enable("betaCat1")
    shinyjs::enable("betaCat2")
    shinyjs::enable("beta_covariates_bin")
    shinyjs::enable("beta_covariatesOptions_bin")
    shinyjs::enable("beta_chooseMethod_bin")
    
    
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  ######################################
  # Taxa Cross-Sectional Data Analysis #
  ######################################
  observeEvent(input$taxa_runbtn_bin,{
    
    validate(
      if (input$covariates_taxa == "Covariate(s)" & is.null(input$covariatesOptions_taxa)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("taxa_runbtn_bin")
    shinyjs::disable("dataType_taxa")
    shinyjs::disable("primvar_taxa")
    shinyjs::disable("taxaCat1")  
    shinyjs::disable("taxaCat2")
    shinyjs::disable("covariates_taxa")
    shinyjs::disable("covariatesOptions_taxa")
    shinyjs::disable("include_species.dend")
    shinyjs::disable("chooseMethod_taxa")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        taxa.categors = c(taxa.categos$cat1, taxa.categos$cat2)
        taxa.bin_categos = c(input$taxaCat1, input$taxaCat2)
        
        rename.cats_ref = taxa.bin_categos[which(taxa.categors == taxa.categos$cat1)]
        rename.cats_com = taxa.bin_categos[which(taxa.categors != taxa.categos$cat1)]
        
        taxa.bin.cat.ref.ori.out <- taxa.bin.cat.ref.ori.func(chooseData$sam.dat, input$primvar_taxa)
        taxa.results$lib.size <- taxa.results$lib.size[names(taxa.results$lib.size) %in% rownames(chooseData$sam.dat)]
        
        if (input$include_species.dend == "Phylum - Species") {
          include = TRUE
        } else {
          include = FALSE
        }
        
        incProgress(1/10, message = "Examining Data in progress")
        sam_dat <- taxa.bin.cat.recode.func(chooseData$sam.dat, input$primvar_taxa, taxa.bin.cat.ref.ori.out,
                                            rename.cats_ref, rename.cats_com)
        
        if (input$covariates_taxa == "None") {
          if (input$chooseMethod_taxa == "Negative binomial regression") {
            taxa.bin.out <- taxa.bin.cat.ref.united.func(input$primvar_taxa, rename.cats_ref,
                                                         rename.cats_com, sam_dat, taxa = chooseData$taxa.out[["count"]])
            
          } else {
            taxa.bin.out <- taxa.bin.cat.ref.united.func(input$primvar_taxa, rename.cats_ref,
                                                         rename.cats_com, sam_dat, taxa = chooseData$taxa.out[[taxa.types$dataType]])
          }
          taxa.results$bin.var <- taxa.bin.out$bin.var
          taxa.results$taxa <- taxa.bin.out$taxa
          
        } else if (input$covariates_taxa == "Covariate(s)") {
          if (input$chooseMethod_taxa == "Negative binomial regression") {
            
            taxa.bin.cov.out <- taxa.bin.cov.cat.ref.united.func(input$primvar_taxa, rename.cats_ref, 
                                                                 rename.cats_com, input$covariatesOptions_taxa, 
                                                                 sam_dat, chooseData$taxa.out[["count"]])
            taxa.bin.cov.out <- taxa_no_miss_cov(sam_dat, taxa.bin.cov.out, input$covariatesOptions_taxa)
          } else {
            taxa.bin.cov.out <- taxa.bin.cov.cat.ref.united.func(input$primvar_taxa, rename.cats_ref, 
                                                                 rename.cats_com, input$covariatesOptions_taxa, 
                                                                 sam_dat, chooseData$taxa.out[[taxa.types$dataType]])
            
            taxa.bin.cov.out <- taxa_no_miss_cov(sam_dat, taxa.bin.cov.out, input$covariatesOptions_taxa)
          }
          taxa.results$bin.var <- taxa.bin.cov.out$bin.var
          taxa.results$taxa <- taxa.bin.cov.out$taxa
          taxa.results$cov.var <- taxa.bin.cov.out$cov.var
          
        }
        taxa.results$taxa.bin.sum.out <- taxa.bin.sum.united.func(taxa.results$bin.var, taxa.results$taxa)
        
        taxa_dataBinvar <- taxa.results$bin.var
        taxa_dataTaxa <- taxa.results$taxa
        
        if (input$chooseMethod_taxa == "Welch t-test" | input$chooseMethod_taxa == "Wilcoxon rank-sum test (default)") {
          if (input$chooseMethod_taxa == "Welch t-test") {
            incProgress(5/10, message = "Welch t-test")
            taxa.t.test.out <- taxa.bin.t.test.united(taxa.results$bin.var, taxa.results$taxa)
            taxa.t.test.q.out <- bin.q.united.func(taxa.t.test.out, method = "BH")
            
            taxa.outputs$DAoutput = taxa.t.test.q.out
            
            nrow = numeric()
            for (r in 1:6) {
              row.num <- ceiling(sum(taxa.t.test.q.out[[r]]$Q.value < 0.05)/4)
              if (row.num > 0) {
                nrow[r] <- row.num
              } else {
                nrow[r] <- 1
              } 
            }
          } else if (input$chooseMethod_taxa == "Wilcoxon rank-sum test (default)") {
            incProgress(5/10, message = "Wilcoxon rank-sum test")
            taxa.wilcox.test.out <- taxa.bin.wilcox.test.united(taxa.results$bin.var, taxa.results$taxa)
            taxa.wilcox.test.q.out <- bin.q.united.func(taxa.wilcox.test.out, method = "BH")
            taxa.wilcox.test.est.added <- taxa.wilcox.test.est.func(taxa.results$bin.var, taxa.results$taxa, rename.cats_ref, rename.cats_com, taxa.wilcox.test.q.out)
            
            taxa.outputs$DAoutput = taxa.wilcox.test.est.added
            
            nrow = numeric()
            for (r in 1:6) {
              row.num <- ceiling(sum(taxa.outputs$DAoutput[[r]]$Q.value < 0.05)/4)
              if (row.num > 0) {
                nrow[r] <- row.num
              } else {
                nrow[r] <- 1
              }
            }
          }
          
          
          incProgress(3/10, message = "Displaying Results in progress")
          output$taxa_display_results = renderUI({
            tagList(
              tabBox(title = strong("Box Plot", style = "color:black"), width = NULL,
                     tabPanel("Phylum", align = "center",
                              plotOutput("rank1", height = nrow[1]*250, width = 750),
                     )
                     ,
                     tabPanel("Class", align = "center",
                              plotOutput("rank2", height = nrow[2]*250, width = 750),
                     )
                     ,tabPanel("Order", align = "center",
                               plotOutput("rank3", height = nrow[3]*250, width = 750),
                     )
                     ,tabPanel("Family", align = "center",
                               plotOutput("rank4", height = nrow[4]*250, width = 750),
                     )
                     ,tabPanel("Genus", align = "center",
                               plotOutput("rank5", height = nrow[5]*250, width = 750),
                     )
                     ,tabPanel("Species", align = "center",
                               plotOutput("rank6", height = nrow[6]*250, width = 750),
                     )
              )
            )
          })
          
          output$rank1 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 1, TRUE)  
          })
          
          output$rank2 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 2, TRUE)  ####all.t.test.united should be used
          })
          
          output$rank3 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 3, TRUE)  ####all.t.test.united should be used
          })
          
          output$rank4 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 4, TRUE)  ####all.t.test.united should be used
          })
          
          output$rank5 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 5, TRUE)  ####all.t.test.united should be used
          })
          
          output$rank6 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 6, TRUE)  ####all.t.test.united should be used
          })
          
          
          taxa.sig <- taxa.sig.dend(taxa.outputs$DAoutput, chooseData$NAadded$tax.tab, "twopi", include)
          
          
          flow.text <- taxa.sig$flow.text
          taxon.tab <- taxa.sig$taxon.tab
          ci.tab.all <- taxa.sig$ci.tab.all
          
          if ( include == FALSE){
            taxon.tab <- taxon.tab[ !grepl("s_", taxon.tab$Taxon, fixed = TRUE), ]
          }
          
          if ( length(ci.tab.all) > 1 ){
            
            for( i in 1:nrow(taxon.tab)){
              if ( ci.tab.all[-1][i] < 0){
                taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "blue")
                taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
              }
              else{
                taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "red")
                taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
              }
              
            }
          }
          
          N <- dim(taxon.tab)[1]
          itr <- ceiling(N/5)
          tab.one   <- data.frame( matrix(ncol=2,nrow=0) )
          tab.two   <- data.frame( matrix(ncol=2,nrow=0) )
          tab.three <- data.frame( matrix(ncol=2,nrow=0) )
          tab.four  <- data.frame( matrix(ncol=2,nrow=0) )
          tab.five  <- data.frame( matrix(ncol=2,nrow=0) )
          
          colnames(tab.one) <- c("ID", "Taxon")
          colnames(tab.two) <- c("ID", "Taxon")
          colnames(tab.three) <- c("ID", "Taxon")
          colnames(tab.four) <- c("ID", "Taxon")
          colnames(tab.five) <- c("ID", "Taxon")
          
          if ( dim(taxon.tab)[1] > 0 ) {
            
            i = 1
            j = 1
            N.rmd <- N %% 5
            N.fix <- N + 5 - N.rmd
            while ( i <= itr) {
              tab.one  [i,] <- taxon.tab[j,]
              tab.two  [i,] <- taxon.tab[j+1,]
              tab.three[i,] <- taxon.tab[j+2,]
              tab.four [i,] <- taxon.tab[j+3,]
              tab.five [i,] <- taxon.tab[j+4,]
              i <- i + 1  
              j <- j + 5
            }
            row.names(tab.one) <- NULL
            row.names(tab.two) <- NULL
            row.names(tab.three) <- NULL
            row.names(tab.four) <- NULL
            row.names(tab.five) <- NULL
            
            tab.one <- na.omit(tab.one)
            tab.two <- na.omit(tab.two)
            tab.three <- na.omit(tab.three)
            tab.four <- na.omit(tab.four)
            tab.five <- na.omit(tab.five)
          }
          
          
          output$taxa_display_dend = renderUI({
            
            box(title = strong("Dendrogram"), width = 12, status = "primary", solidHeader = TRUE,
                
                fluidRow(width = 12, align = "center",
                         div(style = "display: inline-block:vertical-align:top;", grVizOutput("dendrogram", height = 1000, width = 1000)) ),
                br(),
                fluidRow(width = 12, align = "center",
                         tagList(
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("M2sig_1_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("M2sig_2_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("M2sig_3_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("M2sig_4_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("M2sig_5_taxonlist") )
                         )
                )
            )
          })
          
          output$dendrogram = renderGrViz({
            flow.text
          })
          output$M2sig_1_taxonlist <- renderText({
            sig.tab1 <- kable(tab.one, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab1
          })
          output$M2sig_2_taxonlist <- renderText({
            
            sig.tab2 <- kable(tab.two, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab2
          })
          output$M2sig_3_taxonlist <- renderText({
            
            sig.tab3 <- kable(tab.three, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab3
          })
          output$M2sig_4_taxonlist <- renderText({
            
            sig.tab4 <- kable(tab.four, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab4
          })
          output$M2sig_5_taxonlist <- renderText({
            
            sig.tab5 <- kable(tab.five, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab5
          })
          
          
          # 
          # output$M2sig_taxonlist <- renderText({
          #   #sig.taxon.tab
          #   kable(taxon.tab, 'html', booktabs =TRUE, escape = FALSE) %>%
          #     kable_styling(latex_options = c('hold_position')) %>% 
          #     scroll_box(width = "338px", height = "1000px")
          # })
        } 
        
        else if (input$chooseMethod_taxa == "Linear regression"| input$chooseMethod_taxa == "Negative binomial regression" | input$chooseMethod_taxa == "Beta regression") {
          if (input$covariates_taxa == "None") {
            if (input$chooseMethod_taxa =="Linear regression") {
              incProgress(5/10, message = "Linear regression")
              taxa.lm.bin.out <- taxa.bin.lm.united.func(bin.var = taxa.results$bin.var, taxa = taxa.results$taxa)
              taxa.lm.bin.q.out <- bin.q.united.func(taxa.lm.bin.out, method = "BH")
              
              taxa.outputs$DAoutput = taxa.lm.bin.q.out
              nrow <- taxa.forest.plot.pages(taxa.lm.bin.q.out, species.include = include)
              
            } else if (input$chooseMethod_taxa == "Negative binomial regression") {
              incProgress(5/10, message = "Negative binomial regression")
              taxa.bin.glm.nb.q.out <- all.taxa.bin.glm.nb(taxa.results$bin.var, taxa.results$taxa, taxa.results$lib.size)
              
              taxa.outputs$DAoutput <- taxa.bin.glm.nb.q.out
              nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
              
            } else if (input$chooseMethod_taxa == "Beta regression") {
              incProgress(5/10, message = "Beta regression")
              taxa.bin.beta.q.out <- all.taxa.bin.beta(taxa.results$bin.var, taxa.results$taxa)
              
              taxa.outputs$DAoutput <- taxa.bin.beta.q.out
              nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
            }
            
          } else if (input$covariates_taxa == "Covariate(s)") {
            if (input$chooseMethod_taxa =="Linear regression") {
              incProgress(5/10, message = "Linear regression")
              
              taxa.lm.bin.cov.out <- taxa.bin.cov.lm.united.func(bin.var = taxa.results$bin.var, 
                                                                 cov.var = taxa.results$cov.var,
                                                                 taxa = taxa.results$taxa)
              
              taxa.lm.bin.cov.q.out <- bin.q.united.func(taxa.lm.bin.cov.out, method = "BH")
              
              taxa.outputs$DAoutput = taxa.lm.bin.cov.q.out
              nrow <- taxa.forest.plot.pages(taxa.lm.bin.cov.q.out, species.include = include)
              
            } else if (input$chooseMethod_taxa == "Negative binomial regression") {
              incProgress(5/10, message = "Negative binomial regression")
              taxa.bin.cov.glm.nb.q.out <- all.taxa.bin.cov.glm.nb(taxa.results$bin.var, taxa.results$cov.var, taxa.results$taxa, taxa.results$lib.size)
              
              taxa.outputs$DAoutput <- taxa.bin.cov.glm.nb.q.out
              nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
              
            } else if (input$chooseMethod_taxa == "Beta regression") {
              incProgress(5/10, message = "Beta regression")
              taxa.bin.cov.beta.q.out <- all.taxa.bin.cov.beta(taxa.results$bin.var, taxa.results$cov.var, taxa.results$taxa)
              
              taxa.outputs$DAoutput <- taxa.bin.cov.beta.q.out
              nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
            }
          }
          
          forestplot.data <- taxa.forest.plot.pages1(taxa.outputs$DAoutput, chooseData$taxa.names.out, report.type = "Est", species.include = include)
          
          if (any(!is.na(unlist(chooseData$taxa.names.out$duplicates)))) {
            duplicate.taxa <- sapply(strsplit(unlist(chooseData$taxa.names.out$duplicates), " :"),  "[", 1)
            taxon.inplot <- unlist(lapply(forestplot.data$all.text.tab, `[`, i =, j = 3))
            duplicate.texts <- sum(duplicate.taxa %in% taxon.inplot)
          } else {
            duplicate.texts <- 0
          }
          
          if (duplicate.texts>0) {
            output$taxa_display_results = renderUI({
              tagList(
                do.call(tabsetPanel, lapply(1:nrow, function(i) {
                  tabPanel(title = paste0("Page ", i), align = "center",
                           plotOutput(paste0("M2forest", i), height = 800, width = 750),
                           plotOutput(paste0("duplicates", i), height = 11.5*duplicate.texts+10, width = 750))
                })) 
              )
            })
          } else {
            output$taxa_display_results = renderUI({
              tagList(
                do.call(tabsetPanel, lapply(1:nrow, function(i) {
                  tabPanel(title = paste0("Page ", i), align = "center",
                           plotOutput(paste0("M2forest", i), height = 800, width = 750))
                })) 
              )
            })
          }
          
          lapply(1:nrow, function(j) {
            output[[paste0("M2forest", j)]] <- renderPlot({
              taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data, page = j)  
            })
          })
          
          lapply(1:nrow, function(j) {
            output[[paste0("duplicates", j)]] <- renderPlot({
              duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
            })
          })
          
          taxa.sig <- taxa.sig.dend(taxa.outputs$DAoutput, chooseData$NAadded$tax.tab, "twopi", include)
          
          flow.text <- taxa.sig$flow.text
          taxon.tab <- taxa.sig$taxon.tab
          ci.tab.all <- taxa.sig$ci.tab.all
          
          
          if ( include == FALSE){
            taxon.tab <- taxon.tab[ !grepl("s_", taxon.tab$Taxon, fixed = TRUE), ]
          }
          
          if ( length(ci.tab.all) > 1 ){
            
            for( i in 1:nrow(taxon.tab)){
              if ( ci.tab.all[-1][i] < 0){
                taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "blue")
                taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
              }
              else{
                taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "red")
                taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
              }
              
            }
          }
          
          N <- dim(taxon.tab)[1]
          itr <- ceiling(N/5)
          tab.one   <- data.frame( matrix(ncol=2,nrow=0) )
          tab.two   <- data.frame( matrix(ncol=2,nrow=0) )
          tab.three <- data.frame( matrix(ncol=2,nrow=0) )
          tab.four  <- data.frame( matrix(ncol=2,nrow=0) )
          tab.five  <- data.frame( matrix(ncol=2,nrow=0) )
          
          colnames(tab.one) <- c("ID", "Taxon")
          colnames(tab.two) <- c("ID", "Taxon")
          colnames(tab.three) <- c("ID", "Taxon")
          colnames(tab.four) <- c("ID", "Taxon")
          colnames(tab.five) <- c("ID", "Taxon")
          
          if ( dim(taxon.tab)[1] > 0 ) {
            i = 1
            j = 1
            N.rmd <- N %% 5
            N.fix <- N + 5 - N.rmd
            while ( i <= itr) {
              tab.one  [i,] <- taxon.tab[j,]
              tab.two  [i,] <- taxon.tab[j+1,]
              tab.three[i,] <- taxon.tab[j+2,]
              tab.four [i,] <- taxon.tab[j+3,]
              tab.five [i,] <- taxon.tab[j+4,]
              i <- i + 1  
              j <- j + 5
            }
            row.names(tab.one) <- NULL
            row.names(tab.two) <- NULL
            row.names(tab.three) <- NULL
            row.names(tab.four) <- NULL
            row.names(tab.five) <- NULL
            
            tab.one <- na.omit(tab.one)
            tab.two <- na.omit(tab.two)
            tab.three <- na.omit(tab.three)
            tab.four <- na.omit(tab.four)
            tab.five <- na.omit(tab.five)
          }
          
          output$taxa_display_dend = renderUI({
            
            box(title = strong("Dendrogram"), width = 12, status = "primary", solidHeader = TRUE,
                
                fluidRow(width = 12, align = "center",
                         div(style = "display: inline-block:vertical-align:top;", grVizOutput("dendrogram", height = 1000, width = 1000)) ),
                br(),
                fluidRow(width = 12, align = "center",
                         tagList(
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("M2sig_1_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("M2sig_2_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("M2sig_3_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("M2sig_4_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("M2sig_5_taxonlist") )
                         )
                )
            )
          })
          
          output$dendrogram = renderGrViz({
            flow.text
          })
          
          output$M2sig_1_taxonlist <- renderText({
            sig.tab1 <- kable(tab.one, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab1
          })
          
          output$M2sig_2_taxonlist <- renderText({
            sig.tab2 <- kable(tab.two, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab2
          })
          
          output$M2sig_3_taxonlist <- renderText({
            sig.tab3 <- kable(tab.three, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab3
          })
          
          output$M2sig_4_taxonlist <- renderText({
            sig.tab4 <- kable(tab.four, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab4
          })
          
          output$M2sig_5_taxonlist <- renderText({
            sig.tab5 <- kable(tab.five, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab5
          })
          
          output$duplicates = renderPlot({
            duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
          })
          
        }
        
        incProgress(1/10, message = "Displaying Results in progress")
        output$downloadTable_taxa = renderUI({
          tagList(
            p(" ", style = "margin-top: 20px;"),
            box(title = strong("Download Output Table"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"),
                h5("Summary Statistics"),
                downloadButton("tdownloadTabl1", "Download", width = '50%', style = "background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("tdownloadTabl2", "Download", width = '50%', style = "background-color: red3"),
                h5("Dendrogram"),
                downloadButton("gdownload", "Download", width = '50%', style = "background-color: red3")
            )
          )
        })
        
        output$tdownloadTabl1 <- downloadHandler(
          filename = function() {
            paste("Taxa.Sum.Table.zip")
          },
          content = function(sum.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=sum.file, files=dataFiles)
          }
        )
        # output$tdownloadTabl1 <- downloadHandler(
        #   filename = function() {
        #     paste("Taxa.Sum.Table", ".csv", sep="")
        #   },
        #   content = function(file) {
        #     write.list(taxa.results$taxa.bin.sum.out, file, row.names = TRUE)
        #     
        #   }
        # )
        
        output$tdownloadTabl2 <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Output.zip")
          },
          content = function(DA.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.outputs$DAoutput$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=DA.file, files=dataFiles)
          }
        )
        
        # output$tdownloadTabl2 <- downloadHandler(
        #   filename = function() {
        #     paste("Taxa.Analysis.Output", ".csv", sep="")
        #   },
        #   content = function(file) {
        #     write.list(taxa.outputs$DAoutput, file, row.names = TRUE)
        #   }
        # )
        
        output$gdownload <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Graphical.Output", ".html", sep="")
          },
          content = function(file) {
            htmlwidgets::saveWidget(as_widget(taxa.sig.dend(taxa.outputs$DAoutput, chooseData$NAadded$tax.tab, "twopi", include)), file)
          }
        )
        
        ref_string = REFERENCE_CHECK(data_transform = input$dataType_taxa, method_name = isolate(input$chooseMethod_taxa), FDR = isolate("Yes"))
        if (is.null(ref_string)) {
          shinyjs::hide("referencesM2.taxa")
        } else {
          shinyjs::show("referencesM2.taxa")
          output$referencesM2.taxa = renderUI({
            tagList(
              box(title = strong("References"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        shinyjs::enable("taxa_runbtn_bin")
        shinyjs::enable("dataType_taxa")
        shinyjs::enable("primvar_taxa")
        shinyjs::enable("taxaCat1")  
        shinyjs::enable("taxaCat2")
        shinyjs::enable("covariates_taxa")
        shinyjs::enable("covariatesOptions_taxa")
        shinyjs::enable("include_species.dend")
        shinyjs::enable("chooseMethod_taxa")
        
      })
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
}

shinyApp(ui = ui, server = server)
