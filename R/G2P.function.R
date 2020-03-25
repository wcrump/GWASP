#' Simulate a phenotype
#'
#' @param X A genotype matrix (rows = taxa, columns = markers).
#' @param h2 The heritability of the simulated trait.
#' @param alpha The p parameter in an approximated geometric distribution
#' @NQTN The number markers contributing to the phenotype
#' @a2 The proportion of non-additive interactive variance (FIXED at 0 for our current uses)
#' @seed A seed number to fix the randomization (if desired)
#' @return A list of the marker effects (1*NQTN), the phenotype (n*1), the combined additive effects (n*1), the residual effects (n*1), the additive marker indices (1*m), and the interaction marker positions (1*nint)

G2P <- function(X = NULL,h2 = 0.5,alpha = 1,NQTN = 10,distribution = "normal",a2=0){

	n=nrow(X) # number of taxa/samples
	m=ncol(X) # number of markers
	#Sampling QTN
	QTN.position=sample(m,NQTN,replace=F) # Randomly sample markers
	SNPQ=as.matrix(X[,QTN.position]) # Extract the selected markers
	QTN.position
	#QTN effects
	if(distribution=="normal")
	{addeffect=rnorm(NQTN,0,1)
	}else
	{addeffect=alpha^(1:NQTN)} # Simulate with the pseudo-geometric distribution (THIS COULD BE ANY DISTRIBUTION YOU WANT)
	#Simulate phenotype
	effect=SNPQ%*%addeffect # Use matrix multiplication to get the total genetic effect
	effectvar=var(effect)

	#Interaction (DOESN'T appear to be used....)
	cp=0*effect
	nint=4 # The number of interactions appears to be fixed at 4. MAGIC NUMBER???
	if(a2>0&NQTN>=nint){
		for(i in nint:nint){
			Int.position=sample(NQTN,i,replace=F)
			cp=apply(SNPQ[,Int.position],1,prod)
		}
		print(dim(cp))
		cpvar=var(cp)
		intvar=(effectvar-a2*effectvar)/a2
		if(is.na(cp[1]))stop("something wrong in simulating interaction")

		if(cpvar!=0){
			print(c(effectvar,intvar,cpvar,var(cp),a2))
			print(dim(cp))
			cp=cp/sqrt(cpvar)
			cp=cp*sqrt(intvar)
			effectvar=effectvar+intvar
		}else{cp=0*effect}
	}

	residualvar=(effectvar-h2*effectvar)/h2 # Calculate the variance of the residual (scaled to the already simulated effects)
	residual=rnorm(n,0,sqrt(residualvar)) # Sample residual from a normal distribution
	y=effect+residual+cp # Combine simulated effects, residual, and interations
	return(list(addeffect = addeffect, y=y, add = effect, residual = residual, QTN.position=QTN.position, SNPQ=SNPQ,int=cp))
}
