
#######how to install powsimR in R#########
####first download devtools#######
#open terminal ####
####sudo apt-get install r-cran-devtools (terminal cammands)
###then open rstudio###

install.packages("rlang")
install.packages("devtools")
library(devtools)
 ###again open terminal install git##
####sudo apt-get install git-all (install github)####


###In rstudio#####
devtools::install_github("bvieth/powsimR", build_vignettes = TRUE, dependencies = FALSE)
devtools::install_github("bvieth/powsimR")
library("powsimR")

###then run the three different packages installation steps#######

ipak <- function(pkg, repository = c("CRAN", "Bioconductor", "github")) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    # new.pkg <- pkg
    if (length(new.pkg)) {
        if (repository == "CRAN") {
            install.packages(new.pkg, dependencies = TRUE)
        }
        if (repository == "Bioconductor") {
            if (strsplit(version[["version.string"]], " ")[[1]][3] > "4.0.0") {
                if (!requireNamespace("BiocManager")) {
                  install.packages("BiocManager")
                }
                BiocManager::install(new.pkg, dependencies = TRUE, ask = FALSE)
            }
            if (strsplit(version[["version.string"]], " ")[[1]][3] < "3.6.0") {
                stop(message("powsimR depends on packages and functions that are only available in R 4.0.0 and higher."))
            }
        }
        if (repository == "github") {
            devtools::install_github(new.pkg, build_vignettes = FALSE, force = FALSE,
                dependencies = TRUE)
        }
    }
}

####### CRAN PACKAGES#########

cranpackages <- c("broom", "cobs", "cowplot", "data.table", "doParallel", "dplyr",
    "DrImpute", "fastICA", "fitdistrplus", "foreach", "future", "gamlss.dist", "ggplot2",
    "ggpubr", "ggstance", "grDevices", "grid", "Hmisc", "kernlab", "MASS", "magrittr",
    "MBESS", "Matrix", "matrixStats", "mclust", "methods", "minpack.lm", "moments",
    "msir", "NBPSeq", "nonnest2", "parallel", "penalized", "plyr", "pscl", "reshape2",
    "Rmagic", "rsvd", "Rtsne", "scales", "Seurat", "snow", "sctransform", "stats",
    "tibble", "tidyr", "truncnorm", "VGAM", "ZIM", "zoo")



########### BIOCONDUCTOR#####

biocpackages <- c("bayNorm", "baySeq", "BiocGenerics", "BiocParallel", "DESeq2",
    "EBSeq", "edgeR", "IHW", "iCOBRA", "limma", "Linnorm", "MAST", "monocle", "NOISeq",
    "qvalue", "ROTS", "RUVSeq", "S4Vectors", "scater", "scDD", "scde", "scone", "scran",
    "SCnorm", "SingleCellExperiment", "SummarizedExperiment", "zinbwave")


# GITHUB########

githubpackages <- c("cz-ye/DECENT", "nghiavtr/BPSC", "mohuangx/SAVER", "statOmics/zingeR",
    "Vivianstats/scImpute")



                                                                    ##########START################
###power analysis with powsimR#####

####1 Estimate simulation parameter value##### 

arab_data <- readRDS(file.path(getwd(), "Downloads", "R_data", "ch1", "arabidopsis.RDS" ))
means_mock <- rowMeans(arab_data[, c("mock1", "mock2", "mock3")])
means_hrcc <- rowMeans(arab_data[, c("hrcc1", "hrcc2", "hrcc3")])
log2fc <- log2(means_hrcc / means_mock)


prop_de <- sum(abs(log2fc) > 2) / length(log2fc) 
prop_de
###2 Examine the distribution of the log2 fold change ratios#####

finite_log2fc <-log2fc[is.finite(log2fc)]
plot(density(finite_log2fc))
install.packages("extRemes")
extRemes::qqnorm(finite_log2fc )

####3 set up parameter value for the simulation run#####
library(powsimR)
library(dplyr)
install.packages("params")

  params <- estimateParam(
    countData = arab_data,
    Distribution = "NB",
    RNAseq = "bulk",
    normalisation = "TMM" # edgeR method, can be others
  )
  
  de_opts <- DESetup(ngenes=1000,
                     nsims=25,
                     p.DE = prop_de,
                     pLFC= finite_log2fc,
                     sim.seed = 58673
  )

  sim_opts <- SimSetup(
    desetup = de_opts,
    params = params
  )
  
  num_replicates <- c(2, 3, 5, 8, 12,15)
####4 Run the simulation######
  simDE <- simulateDE(n1 = num_replicates,
                      n2 = num_replicates,
                      sim.settings = sim_opts,
                      DEmethod = "edgeR-LRT",
                      normalisation = "TMM",
                      verbose = FALSE)

####5 Run the evalution of the simulation#########

   evalDE <- evaluateDE(simRes = simDE,
                                 alpha.type = 'adjusted',
                                 MTC = 'BH',
                                 alpha.nominal = 0.1,
                                 stratify.by = 'mean',
                                 filter.by = 'none',
                                 strata.filtered = 1,
                                 target.by = 'lfc',
                                 delta = 0)

####6 Plot the evolution#######
  plotEvalDE(evalRes = evalDE,
                    rate='marginal',
                    quick=FALSE, annot=TRUE)


log2fc_func <- function(x){ rnorm(x, 0, 2)} 
prop_de = 0.1

de_opts <- DESetup(ngenes=1000,
                     nsims=25,
                     p.DE = prop_de,
                     pLFC= log2fc_func,
                     sim.seed = 58673
  )
