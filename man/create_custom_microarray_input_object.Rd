% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_input.R
\name{create_custom_microarray_input_object}
\alias{create_custom_microarray_input_object}
\title{Create a generic microarray based input object}
\usage{
create_custom_microarray_input_object(diff_genes = NULL,
  limma_probe_table = NULL, annotated_exprs_matrix = NULL,
  expression_matrix = NULL, annotation_table = NULL,
  group1_indici = NULL, group2_indici = NULL, group1_label = NULL,
  group2_label = NULL, settings = NULL)
}
\arguments{
\item{diff_genes}{A 2 two column data.frame where the first column are genes and the second column are p-values}

\item{limma_probe_table}{A data.frame from \code{limma topTable} with added gene annotation}

\item{annotated_exprs_matrix}{A matrix where the rows are genes and the columns samples.}

\item{expression_matrix}{Normalized expression matrix where the samples are columns and probes are rows}

\item{annotation_table}{A dataframe providing annotation for the probes. The dataframe should have 3 columns: 
\itemize{
\item {PROBEID}: The probe id as it is in the expression matrix
\item{IDENTIFIER}: The entrez id (if available) associated with the probe
}}

\item{group1_indici}{vector containing indici for samples belonging to group 1 (Column numbers)}

\item{group2_indici}{vector containing indici for samples belonging to group 2 (Column numbers)}

\item{group1_label}{Label for each group 1, for example "patient" or "control"}

\item{group2_label}{Label for each group 2, for example "patient" or "control"}

\item{settings}{Settings used to generate the object. Used only internally by the package.}
}
\value{
The function returns an object of class "MODifieR_input". The object is a named list containing the
following components:
\item{diff_genes}{A 2 two column data.frame where the first column are genes and the second column p-values}
\item{limma_probe_table}{A data.frame from \code{limma topTable} with added gene annotation}
\item{annotated_exprs_matrix}{A matrix where the rows are genes and the columns samples. Probes have been collapsed
into genes using \code{collapse_method}}
\item{expression_matrix}{A matrix, the original input expression matrix}
\item{annotation_table}{A data.frame, the original annotation table used to annotate the probes}
\item{group_indici}{A named list containing 2 numeric vectors. The names are the group labels and the values 
are the group indici}
}
\description{
Create a generic microarray based input object
}
\details{
This function allows the creation of a generic microarray based input object with the same class as objects created by
\code{\link{create_input_microarray}}. This can be useful in the cases where you already have differentially 
expressed genes, an annotated expression matrix or both and want to wrap that into an input object
to use in downstream analysis. All arguments are optional.
}
\seealso{
\code{\link{create_input}}
}
\author{
Dirk de Weerd
}
