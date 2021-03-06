% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manhattan.plot.function.R
\name{manhattan_plot}
\alias{manhattan_plot}
\title{Generate a Manhattan Plot of a GWAS result}
\usage{
manhattan_plot(marker_map, pvals, QTN_index = c(), trait = "unknown")
}
\arguments{
\item{marker_map}{A data.frame containing marker 'rs', 'chr', and 'pos' (should be in the same order as the GWAS results)}

\item{pvals}{The untransformed p-values from a GWAS test}

\item{QTN_index}{The indices of any known QTNs}

\item{trait}{The name of the trait}
}
\value{
A ggplot object Manhattan plot
}
\description{
Generate a Manhattan Plot of a GWAS result
}
