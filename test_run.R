#bring in the data
gen.data <- read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
gen.map <- read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
phen.data <- read.table(file="http://zzlab.net/GAPIT/data/CROP545_Phenotype.txt",head=T)
covs <- read.table(file="http://zzlab.net/GAPIT/data/CROP545_Covariates.txt",head=T)

library(devtools)
install_github("wcrump/GWASP", ref = "master",
			auth_token = github_pat(quiet), host = "https://api.github.com", quiet = FALSE, force = T)
library(GWASP)

pvals <- GLM.func(geno = gen.data, pheno = phen.data, covariates = covs, PCs = 3, thresh = 0.2)
