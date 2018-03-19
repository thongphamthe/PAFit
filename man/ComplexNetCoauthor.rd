\name{Coauthorship network of scientists in the field of complex networks}
\docType{data}
\alias{ComplexNetCoauthor}
\alias{coauthor.net}
\alias{coauthor.truetime}
\alias{coauthor.author_id}
\title{A collaboration network between authors of papers in the field of complex networks with article time-stamps}
\description{
The dataset is collaboration network of authors of network science articles with article time-stamps. An edge between two authors represents an article in common. Time stamps denote article publication dates. The network without time-stamps was compiled by Mark Newman in May 2006 from the bibliographies of two review articles on networks, M. E. J. Newman, SIAM Review 45, 167-256 (2003) and S. Boccaletti et al., Physics Reports 424, 175-308 (2006), with a few additional references added by hand. Paul Sheridan independently supplemented the network with time-stamps and some basic metadata in June 2015. The network is undirected with monthly resolution, and contains no duplicated edges. \code{coauthor.net} contains the network. \code{coauthor.truetime} contains the real times of processed time-stamps. Finally \code{coauthor.author_id} contains author names.  

Reference: M. E. J. Newman, Finding community structure in networks using the eigenvectors of matrices, Preprint physics/0605087 (2006).

The original network data that Mark Newman compiled is available for download at \url{http://www-personal.umich.edu/~mejn/netdata/netscience.zip}

The dataset with article time-stamps is available for download at \url{http://www.paulsheridan.net/files/collabnet.zip}

}
\usage{data(ComplexNetCoauthor)}
\format{\code{coauthor.net} is a matrix with 2849 rows and 3 columns. Each row is an edge with the format (author id 1, author id 2, time_stamp). \code{coauthor.truetime} is a two-column matrix whose each row is (time_stamp, real time). \code{coauthor.author_id} is a two-column matrix whose each row is (author id, author name).}
\source{http://www.paulsheridan.net/files/collabnet.zip}