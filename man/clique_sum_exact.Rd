% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clique_sum_exact.R
\name{clique_sum_exact}
\alias{clique_sum_exact}
\title{Clique Sum Exact}
\usage{
clique_sum_exact(MODifieR_input, db, clique_significance = 0.01,
  deg_cutoff = 0.05, min_clique_size = 5, min_deg_in_clique = 3,
  multiple_cores = T, n_cores = 4, dataset_name = NULL)
}
\arguments{
\item{MODifieR_input}{A MODifieR input object produced by one of the \code{create_input} functions}

\item{db}{A clique database created by \code{\link{build_clique_db}}}

\item{clique_significance}{p-value for cliques to be considered significant}

\item{deg_cutoff}{p-value cutoff for differentialy expressed genes}

\item{min_clique_size}{Minimal size for cliques}

\item{min_deg_in_clique}{Minimum number of DEGs to be present in a clique}

\item{multiple_cores}{Parallel process using multiple cores?}

\item{n_cores}{Number of cores to use}

\item{dataset_name}{Optional name for the input object that will be stored in the settings object.
Default is the variable name of the input object}
}
\value{
clique_sum_exact returns an object of class "MODifieR_module" with subclass "Clique_Sum_exact". 
This object is a named list containing the following components:
\item{module_genes}{A character vector containing the genes in the final module}
\item{settings}{A named list containing the parameters used in generating the object}
}
\description{
An implementation of the clique-based disease module inference method proposed by Barrenäs et al.
}
\details{
Clique_sum_exact finds cliques of at least size \code{min_clique_size} that are 
 significantly enriched with DEGs. 
 The union of maximal cliques with a Fisher-exact test p-value below \code{clique_significance} and at 
 least \code{min_deg_in_clique} is the final disease module.
}
\references{
\cite{Barrenäs F, Chavali S, Alves AC, et al. Highly interconnected genes in 
disease-specific networks are enriched for disease-associated polymorphisms. 
Genome Biology. 2012;13(6):R46. doi:10.1186/gb-2012-13-6-r46.}
}
\author{
Dirk de Weerd
}
