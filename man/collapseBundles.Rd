% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collapseBundles.R
\name{collapseBundles}
\alias{collapseBundles}
\title{Collapse bundles of transcripts, discard any that represent pointless tests,
and optionally prune any whose joined bundle IDs tend to choke downstream 
packages for e.g. pathway- or network-based enrichment analysis.  Note that 
this function may or may not be optimal for your RNAseq experiment. Please 
refer to 'Details' for some thought exercises about the nature of 'genes'.}
\usage{
collapseBundles(kexp, bundleID = "gene_id", read.cutoff = 1,
  discardjoined = TRUE)
}
\arguments{
\item{kexp}{A KallistoExperiment (or something very much like it)}

\item{bundleID}{The column (in mcols(features(kexp))) of the bundle IDs}

\item{read.cutoff}{Discard transcripts and bundles with < this many counts}

\item{discardjoined}{Discard bundles with IDs "joined" by a ";"?  (TRUE)}
}
\value{
a matrix of summarized counts per sample bundle
}
\description{
Collapse bundles of transcripts, discard any that represent pointless tests,
and optionally prune any whose joined bundle IDs tend to choke downstream 
packages for e.g. pathway- or network-based enrichment analysis.  Note that 
this function may or may not be optimal for your RNAseq experiment. Please 
refer to 'Details' for some thought exercises about the nature of 'genes'.
}
\details{
This function sums the estimated counts for each transcript within 
a bundle of transcripts (where "bundle" is a user-defined identifier, often 
but not always a 'gene', sometimes a biotype or a class of repeat elements).
The default approach is to discard all rows where the maximum count is less 
than the specified read.cutoff. Since the default cutoff is 1, this means 
discarding transcripts (and bundles) that were not be detected in any sample.
(Filtering tends to increase statistical power at a given false-positive rate
per Bourgon et al, 2010, \link{http://www.pnas.org/content/107/21/9546.long})
}
\seealso{
collapseTranscripts
}

