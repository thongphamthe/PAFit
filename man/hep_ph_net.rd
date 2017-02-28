\name{HEP-PH citation network}
\docType{data}
\alias{hep_ph.rda}
\alias{hep_ph_net}
\title{A citation network between papers on the HEP-PH section of arXiv}
\description{
  The dataset contains a citation network between papers on the HEP-PH section of arXiv. The resolution is yearly.
}
\usage{data(hep_ph)}
\format{The dataset is a matrix with 360159 rows and 3 columns. Each row is an edge with the format (from_paper id, to_paper id, time_stamp). The time-stamps have been converted to labels.}
\source{https://snap.stanford.edu/data/cit-HepPh.html}
\references{
  https://snap.stanford.edu/data/cit-HepPh.html
}