\name{Coauthorship network of scientists in the complex network field}
\docType{data}
\alias{ComplexNetCoauthor.rda}
\alias{coauthor.net}
\alias{coauthor.truetime}
\alias{coauthor.author_id}
\title{A collaboration network between authors of papers in the  complex network fields}
\description{
The dataset contains a coauthorship network of scientists working on network theory and experiment, as compiled by M. Newman in May 2006.  The network was compiled from the bibliographies of two review articles on networks, M. E. J. Newman, SIAM Review 45, 167-256 (2003) and S. Boccaletti et al., Physics Reports 424, 175-308 (2006), with a few additional references added by hand.

The network is undirected with monthly resolution, and contains no duplicated edges. \code{coauthor.net} contains the network. \code{coauthor.truetime} contains the real times of processed time-stamps. Finally \code{coauthor.author_id} contains author names.  

If you make use of these data, please cite M. E. J. Newman, Finding
community structure in networks using the eigenvectors of matrices,
Preprint physics/0605087 (2006).

The original dataset is available for download at www.paulsheridan.net

}
\usage{data(ComplexNetCoauthor)}
\format{\code{coauthor.net} is a matrix with 2849 rows and 3 columns. Each row is an edge with the format (author id 1, author id 2, time_stamp). \code{coauthor.truetime} is a two-column matrix whose each row is (time_stamp, real time). \code{coauthor.author_id} is a two-column matrix whose each row is (author id, author name).}
\source{www.paulsheridan.net}
\references{
  www.paulsheridan.net
}