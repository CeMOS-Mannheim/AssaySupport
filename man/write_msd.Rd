% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_data.R
\name{write_msd}
\alias{write_msd}
\title{write spectra as MSD to disk}
\usage{
write_msd(
  spectra = spectra_align,
  Dir,
  spectraNames = unique(paste0(names(spectra))),
  namePrefix = "",
  peaks = NA,
  annotation_df = NA
)
}
\arguments{
\item{spectra}{List or single \code{MALDIquant::MassSpectrum}, spectra to export.}

\item{Dir}{Character, folder to write spectra}

\item{spectraNames}{Character vector, names of spectra by default the \code{spectra} are expected to be named.}

\item{namePrefix}{Character, prefix to add to filenames}

\item{peaks}{List or single \code{MALDIquant::MassPeaks}, peaks to be added to msd file. Set to NA if no peaks should be added.}

\item{annotation_df}{data.frame containing assigned species. Must contain columns: "species", "mz.theo", "plotIdx", "int", "mz"}
}
\description{
write spectra as MSD to disk
}
