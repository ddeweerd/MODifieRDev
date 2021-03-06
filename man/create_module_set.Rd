% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module_set.R
\name{create_module_set}
\alias{create_module_set}
\title{Consensus module}
\usage{
create_module_set(min_frequency, module_list)
}
\arguments{
\item{min_frequency}{Minimal number of MODifieR modules that a gene should be present in in order
to include it in the final module}

\item{module_list}{A list of MODifieR modules}
}
\value{
Returns an object of class "MODifieR_module" with subclass "Module_set". 
This object is a named list containing the following components:
\item{module_genes}{A character vector containing the genes in the final module}
\item{module_gene_list}{A named list containing the module genes from the original modules}
\item{gene_frequency}{Table containing all the genes present in the modules and their frequency}
\item{method_by_gene}{A named list where the elements are the modules the genes have been found in
and the names are the gene names}
\item{gene_by_method}{A table containing the gene frequencies by combination of methods}
\item{settings}{A named list containing the parameters used in generating the object}
}
\description{
Consensus module
}
\details{
Get a consensus module that is composed of genes present in at least \code{min_frequency} genes.
If the input \code{module_list} is unnamed, the subclasses of the MODifieR objects will be used.
If there is more than 1 of a given subclass present in \code{module_list}, a number will be appended.
}
\author{
Dirk de Weerd
}
