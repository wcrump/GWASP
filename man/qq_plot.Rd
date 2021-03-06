% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qq.plot.function.R
\name{qq_plot}
\alias{qq_plot}
\title{Generate a Quantile-Quantile Plot of a GWAS result}
\usage{
qq_plot(marker_map, pvals, QTN_index, trait = "unknown")
}
\arguments{
\item{marker_map}{A data.frame containing marker 'rs', 'chr', and 'pos' (should be in the same order as the GWAS results)}

\item{pvals}{The untransformed p-values from a GWAS test}

\item{QTN_index}{The indices of any known QTNs}

\item{trait}{The name of the trait}
}
\value{
A ggplot object quantile-quantile plot
}
\description{
Generate a Quantile-Quantile Plot of a GWAS result
}
