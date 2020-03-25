#' Generate a Manhattan Plot of a GWAS result
#'
#' @param marker_map A data.frame containing marker 'rs', 'chr', and 'pos' (should be in the same order as the GWAS results)
#' @param pvals The untransformed p-values from a GWAS test
#' @param QTN_index The indices of any known QTNs
#' @param trait The name of the trait
#' @return A ggplot object Manhattan plot

manhattan_plot <- function(marker_map, pvals, QTN_index = c(), trait = "unknown")
{
	marker_map$pvals <- -log10(pvals) # Add pvalues to the data.frame and log10 transform

	marker_map$comb_pos <- marker_map$chrom * 1e9 + marker_map$pos
	manhattan_plot <- ggplot(marker_map, aes(x = 1:nrow(marker_map), y = pvals, color = factor(chrom))) +
		geom_point() +
		geom_vline(xintercept = QTN_index, color = "red") +
		labs(title = paste("GWAS manhattan plot for trait:", trait),
			y = "-log10(p-value)",
			x = "Marker Position",
			color = "Chromosome")

	return(manhattan_plot)
}
