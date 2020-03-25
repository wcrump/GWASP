## Dependencies ----
# Download the GWASP package from Github: https://github.com/wcrump/GWASP
# Once you have made a copy of the repository your working directory, file paths should work as written
source("R/GLM.function.R") #This is our main function, GWASP
source("R/filter.pca.R") #This function filters PC's calculated in the PCA step to incorporate into our GLM
manhattan_plot <- function(marker_map, pvals, trait = "unknown")
{
	marker_map$pvals <- -log10(pvals) # Add pvalues to the data.frame and log10 transform

	marker_map$comb_pos <- marker_map$chr * 1e9 + marker_map$pos
	manhattan_plot <- ggplot(marker_map, aes(x = 1:nrow(marker_map), y = pvals, color = factor(chr))) +
		geom_point() +
		labs(title = paste("GWAS manhattan plot for trait:", trait),
			y = "-log10(p-value)",
			x = "Marker Position",
			color = "Chromosome")

	return(manhattan_plot)
}
qq_plot <- function(marker_map, pvals, trait = "unknown")
{
	marker_map$pvals <- -log10(pvals) # Add pvalues to the data.frame and log10 transform

	exp_pval_dist <- -log10(runif(nrow(marker_map), 0, 1)) # Sample random p-values from a uniform distribution between 0 and 1

	qq_plot <- ggplot(marker_map, aes(x = sort(exp_pval_dist),
							    y = sort(marker_map$pvals))) +
		geom_point() +
		geom_abline(intercept = 0, slope = 1, color = "red") +
		labs(title = "Q-Q Plot",
			x = "-log10(p-values) Expected",
			y = "-log10(p-values) Observed")

	return(qq_plot)
}
GWASbyCor=function(X,y){
	n=nrow(X)
	r=cor(y,X)
	n=nrow(X)
	t=r/sqrt((1-r^2)/(n-2))
	p=2*(1-pt(abs(t),n-2))
	zeros=p==0
	p[zeros]=1e-10
	return(p)
}
# Download Data
covariate.data <- read.table("Data/CROP545_Covariates.txt", header = T)
phenotype.data <- read.table("Data/CROP545_Phenotype.txt", header = T)
genotype.data <- read.table("Data/mdp_numeric.txt", header = T)
marker.map <- read.table("Data/mdp_SNP_information.txt", header = T)
colnames(marker.map) <- c("rs", "chr", "pos")

## Perform GWAS using GWASP function ----
p.values <- GWASP(X = genotype.data, y = phenotype.data, C = covariate.data, PCs = 3, r = 0.3)

## Create Manhattan and QQ Plots ----
p.vec <- as.vector(p.values) #transforms data type of p-values to run in manhattan plot function
GWASP_manhattan <- manhattan_plot(marker_map = marker.map, pvals = p.vec) #creates manhattan plot

qq_plot(marker_map = marker.map, pvals = p.vec) #creates a QQ plot
#we can see our p-values start to deviate sharply from what we'd expect to see with a normal distribution

## Genome-wide Threshold and Associated SNP's ----
type1=c(0.01, 0.05, 0.1, 0.2) #vector of type 1 error rates corresponding to 1%, 5%, 10%, and 20%
cutoff=quantile(p.vec,type1,na.rm=T) #calculation of quantiles
genome.thresh <- cutoff[1] #with a type 1 error rate of 1%, we can excpect SNPs with p-values higher than this to be associated SNPs

associated.SNP.index <- p.values <= genome.thresh #creates index of associated SNPs
associated.SNPs <- marker.map$pos[associated.SNP.index]

GWASP_manhattan <- manhattan_plot(marker_map = marker.map, pvals = p.vec, QTN_index = associated.SNPs)
#Comments on MAF of associated SNPs

## GWASP is Superior to GWASbyCor ----
