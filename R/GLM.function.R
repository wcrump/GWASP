#' GLM Function
#'
#' @param geno A matrix of genotype data with dimensions n x m
#' @param pheno A matrix of phenotype data with dimensions n x 1
#' @param covariates A matrix of covariate data with dimensions n x t
#' @param PCs A matrix of principle components with dimensions n x ? (variable number of PCs)
#' @param thresh A numeric (0-1) indicating the correlation value of any PC-covariate pair at which the PC will be exluded from model
#' @return A matrix of p values with dimensions 1 x m

GLM.func <- function(geno = NULL, pheno = NULL, covariates = NULL, PCs = 1, thresh = 0.2){
	working.geno <- geno[,-1] #pull user input into genotype matrix, possible data transformation?
	working.pheno <- pheno[,-1] #pull user input into phenotype matrix
	working.cov <- covariates[,-1] #pull user input into covariate matrix
	working.PCs <- PCs #pull in user-defined number of PCs into PC matrix
	working.thresh <- thresh #pull in user-defined correlation threshold between any given PC-covariate pair

	PCA <- prcomp(working.geno) #perform PCA on the genotypes
	filtered.PCs <- filter.pca(PCA = PCA, covs = working.cov, threshold = working.thresh) #filter PCs by covariates according to correlation threshold supplied by user
	used.PCs <- filtered.PCs$x[,1:working.PCs]

	sample.n <- nrow(working.geno) #number of samples in data
	marker.n <- ncol(working.geno) #number of markers in data

	p.matrix <- matrix(data = NA, nrow = 1, ncol = marker.n) #creation of blank matrix with 1 row and marker.n number of columns

	#for loop through genome
	for (i in 1:marker.n) {
		SNP = working.geno[,i] #SNP value set to i column in working.geno
		if(max(SNP)==min(SNP)){ #if monomorphic marker, p = 1
			p.value = 1
		}
		else{
			coeff.matrix <- as.matrix(cbind(1, working.cov, used.PCs, SNP)) #creates matrix of X values to calculate b on
			X.matrix <- t(coeff.matrix)%*%coeff.matrix #creates X^2 multiplication matrix
			inverse.X <- solve(X.matrix) #creates inverse of X^2 matrix
			Y.matrix <- t(coeff.matrix)%*%working.pheno #creates X*Y matrix

			b.effect <- inverse.X%*%Y.matrix #multiple two matrices we just created to calculate b
			yb <- coeff.matrix%*%b.effect #intermediary calculation using matrix multiplication
			e <- working.pheno-yb #calculating residual

			pheno.n <- length(working.pheno)
			var.e <- sum(e^2)/(pheno.n-1) #calculating variance of residuals
			var.t <- inverse.X*var.e #calculating variance of t-stat
			t.stat <- b.effect/sqrt(diag(var.t)) #calculating t-stat using b and var.t

			p.value <- 2*(1-pt(abs(t.stat), pheno.n-2))
		} #monomorphic marker loop
		p.matrix[i] <- p.value[length(p.value)]
	} #SNP loop
	return(p.matrix)
} #end function

