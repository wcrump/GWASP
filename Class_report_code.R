## Dependencies ----
install.packages("ggplot2")
library(ggplot2)
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

## Perform GWAS using GWASP function ----
p.values <- GLM.func(geno = genotype.data, pheno = phenotype.data, covariates = covariate.data, PCs = 3, thresh = 0.3)

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
associated.SNPs <- which(associated.SNP.index) #let's us know which of our 3000 markers are our significant ones

GWASP_manhattan <- manhattan_plot(marker_map = marker.map, pvals = p.vec, QTN_index = associated.SNPs) #create manhattan plot
GWASP_manhattan <- GWASP_manhattan + geom_hline(yintercept = -log10(genome.thresh), linetype = "dashed", color = "black") #add significance threshold
GWASP_manhattan #view manhattan plot

#Comments on MAF of associated SNPs
sig.snps.geno <- genotype.data[,associated.SNPs]

minor_allele_freq <- apply(sig.snps.geno, 2, function(x) #we're applying the function(x) across all columns, designated by the 2, in working_geno_40scaffolds
{
	allele_freq0 <- (sum(x == 0)*2 + sum(x == 1)) / (sum(!is.na(x)) * 2) #frequency of allele 0 is calculated, includes half of heterozygote frequency
	allele_freq2 <- (sum(x == 2)*2 + sum(x == 1)) / (sum(!is.na(x)) * 2) #frequency of allele 2 is calculated, includes half of heterozygote frequency
	return(min(allele_freq0, allele_freq2)) #now we compare allele freq 0 to allele freq 2 and return the lesser of the two
})

ggplot(data = as.data.frame(minor_allele_freq), aes(x= seq_along(minor_allele_freq), y = minor_allele_freq)) + geom_point() +
	labs(title = "Distribution of Minor Allele Frequency of Significant SNPs", x = "SNP Index", y = "Minor Allele Frequency")
## GWASP is Superior to GWASbyCor, not by speed but by statistical power ----
#time an empty for loop to compare our two function times to, repeat 30 times
empty_counter <- 0 #create counter
empty_times <- vector('numeric') #create empty vector to store our loop times in
repeat{
	empty_start <- proc.time() #start timer
	for (i in 1:100){ #empty loop

	}
	empty_end <- proc.time() #end timer
	empty_elapsed <- empty_end[3] - empty_start[3] #calculate elapsed time passed
	empty_times <- append(empty_times, values = empty_elapsed) #append new time to time vector
	empty_counter <- empty_counter +1 #add one to our counter
	if (empty_counter == 30){ #break the loop once our counter reaches 30
		break
	}
}
mean(empty_times) #calculate mean of empty loop
sd(empty_times) #calculate sd of empty loop

#we can use the same code structure to time our function
GWASP_counter <- 0
GWASP_times <- vector('numeric')
repeat{
	GWASP_start <- proc.time()
	p.GWASP <- GLM.func(geno = genotype.data, pheno = phenotype.data, covariates = covariate.data, PCs = 1, thresh = 0.2)
	GWASP_end <- proc.time()
	GWASP_elapsed <- GWASP_end[3] - GWASP_start[3]
	GWASP_times <- append(GWASP_times, values = GWASP_elapsed)
	GWASP_counter <- GWASP_counter +1
	if (GWASP_counter == 30) {
		break
	}
}
mean(GWASP_times)
sd(GWASP_times)

#we want to perform GWASbyCor and see if our GWAS function is faster or slower
GWASbyCor_counter <- 0
GbyCor_times <- vector('numeric')
repeat{
	GbyCor_start <- proc.time()
	p.GbyCor <- GWASbyCor(X = genotype.data[,-1], y = phenotype.data[,-1])
	GbyCor_end <- proc.time()
	GbyCor_elapsed <- GbyCor_end[3] - GbyCor_start[3]
	GbyCor_times <- append(GbyCor_times, values = GbyCor_elapsed)
	GWASbyCor_counter <- GWASbyCor_counter +1
	if (GWASbyCor_counter == 30) {
		break
	}
}
mean(GbyCor_times)
sd(GbyCor_times)
#we can see that our GWAS function is slower than GWASbyCor, so we'll have to beat it on statistical power
#we're going to use the GAPIT function GAPIT.FDR.Type1 to compare power and FDR between the two methods
nrep=30 #number of replicates
set.seed(89760)
GWASP.Rep=replicate(nrep, { #we're using the replicate function to perform GWAS and then calculate power and FDR multiple times
	GWASP.sim <- GLM.func(geno = genotype.data, pheno = phenotype.data, covariates = covariate.data) #perform GWAS
	seqQTN <- associated.SNPs # reminder that we're using our 31 SNPs
	myGWAS <- cbind(marker.map,t(GWASP.sim),NA) #appends our GWAS results to our marker map into a single data frame
	#calculates power and FDR using the GAPIT.FDR.TypeI function, most parameters kept to default
	myStat <- GAPIT.FDR.TypeI(WS=c(1e0,1e3,1e4,1e5), GM = marker.map, seqQTN = associated.SNPs, GWAS = myGWAS, maxOut = 100,MaxBP=1e10)
})

set.seed(89761)
GbyCor.Rep <- replicate(nrep, { #we can use the same replicate function structure to calculate FDR and power for GWASbyCor
	GbyCor.sim <- GWASbyCor(X = genotype.data[,-1], y = phenotype.data[,-1])
	seqQTN <- associated.SNPs
	myGWAS <- cbind(marker.map, t(GbyCor.sim), NA)
	mySTAT <- GAPIT.FDR.TypeI(WS = c(1e0, 1e3, 1e4, 1e5), GM = marker.map, seqQTN = associated.SNPs, GWAS = myGWAS, maxOut = 100, MaxBP = 1e10)
})

power.GWASP <- GWASP.Rep[[2]] #pulling out the power of the function for use in plotting
power.GbyCor <- GbyCor.Rep[[2]] #pulling out the power of the function for use in plotting

#FDR
gwasp.fdr.seq <- seq(3,length(GWASP.Rep),7) #creates index of where FDR values are stored of every replicate
gwasp.fdr <- GWASP.Rep[gwasp.fdr.seq] #pulls all FDR outputs into a new list of just FDR data
gwasp.fdr.mean <- Reduce ("+", gwasp.fdr) / length(gwasp.fdr) #calculates the FDR mean for each individual replicate

gbycor.fdr.seq <- seq(3, length(GbyCor.Rep), 7)
gbycor.fdr <- GbyCor.Rep[gbycor.fdr.seq]
gbycor.fdr.mean <- Reduce ("+", gbycor.fdr) / length(gbycor.fdr)

theColor=rainbow(4) #creates color vector
plot(gwasp.fdr.mean[,1],power.GWASP , type="b", col=theColor [1],xlim=c(0,1)) #plots GWASP FDR means against power
for(i in 2:ncol(gwasp.fdr.mean)){ #shows different 'resolutions' as different colors
	lines(gwasp.fdr.mean[,i], power.GWASP , type="b", col= theColor [i])
}
plot(gbycor.fdr.mean[,1],power.GbyCor , type="b", col=theColor [1],xlim=c(0,1)) #plots GbyCor FDR means against power
for(i in 2:ncol(gbycor.fdr.mean)){ #shows different 'resolutions' as different colors
	lines(gbycor.fdr.mean[,i], power.GbyCor , type="b", col= theColor [i])
}
