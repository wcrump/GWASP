#' Filter PCs
#'
#' Filter out PCs colinear with other covariates.
#' @param PCA The output list from the prcomp() function.
#' @param covs A matrix/data frame of covariates (rows = individuals, columns = covariates)
#' @param threshold A numerical value (0-1) indicating the maximum allowed correlation between any one covariate and any one PC
#' @return A prcomp() list output with PCs filtered

filter.pca <- function(PCA, covs=NULL, threshold=0.2){

	if(is.null(covs)){
		return(PCA)
	}else{
		correlation.matrix <- cor(PCA$x, covs)
		indices <- apply(correlation.matrix, 2, function(x){x>=threshold})
		index.vect <- apply(indices, 2, function(x){which(x)})

		PCA$x <- PCA$x[,-index.vect]

	}
	return(PCA)

}
