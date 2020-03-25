#' Perform GWAS by Correlation
#'
#' @param X A genotype matrix with dimensions n x m (rows = taxa, columns = markers).
#' @param y A phenotype vector of length n
#' @return A vector of p values of length m


#GWAS by correlation
GWASbyCor=function(X,y){
	n=nrow(X)
	r=cor(y,X)
	n=nrow(X)
	t=r/sqrt((1-r^2)/(n-2))
	p=2*(1-pt(abs(t),n-2))
	zeros=p==0
	p[zeros]=1e-10
	return(p)}
