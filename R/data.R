#'Expression matrix
#'
#'An example dataset from GEO (accession number GSE4588) containing systemic lupus erythematosus patients
#'and healthy controls
#'
#' @details 
#'The expression matrix contains 16 microarray samples. Columns 1:9 are healthy controls,
#'columns 10:16 are SLE patients. There are 52307 probes in this dataset.
#'@seealso
#'\code{\link{probe_annotation}}
#'
#'\code{\link{ppi_network}}
#'@references
#'\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4588}.
#'
#'
#'
"expression_matrix"

#'Probe annotation
#'
#'Probe annotation for the example dataset \code{\link{expression_matrix}}.
#'The annotation is taken from Bioconductor package hgu133plus2.db and contains 52307 probes.
#'@seealso
#'\code{\link{expression_matrix}}
#'
#'\code{\link{ppi_network}}
#'@references
#'\url{https://bioconductor.org/packages/release/data/annotation/html/hgu133plus2.db.html}
#'
#' @references \cite{Carlson M (2016). hgu133plus2.db: Affymetrix Human Genome U133 Plus 2.0 
#' Array annotation data (chip hgu133plus2). R package version 3.2.3.} 
#'
"probe_annotation"

#'PPi network
#'
#'A Protein-Protein interaction network from STRING
#'
#'@details
#'This PPi network from STRING is version 7.1 and filtered for interactions with a 
#'confidence score higher than 700. The identifiers for the genes are ENTREZ.
#'There are 64672 interactions (rows) in the dataframe. Each row denotes an interaction;
#'column 1 contains the first gene, column 2 the second. The third column is the confidence score.
#'
#'
#'@seealso
#'\code{\link{probe_annotation}}
#'
#'\code{\link{expression_matrix}}
#'@references
#'\url{https://string-db.org/}
#'
"ppi_network"