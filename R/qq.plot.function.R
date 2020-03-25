#' Generate a Quantile-Quantile Plot of a GWAS result
#'
#' @param marker_map A data.frame containing marker 'rs', 'chr', and 'pos' (should be in the same order as the GWAS results)
#' @param pvals The untransformed p-values from a GWAS test
#' @param QTN_index The indices of any known QTNs
#' @param trait The name of the trait
#' @return A ggplot object quantile-quantile plot

qq_plot <- function(marker_map, pvals, QTN_index, trait = "unknown")
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
