\name{plot.PAFit_net}
\alias{plot.PAFit_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot a \code{PAFit_net} object
}
\description{
  This function plots a \code{PAFit_net} object. There are four options of \code{plot} to specify the type of plot. 
  
  The first two concern plotting the graph in \code{$graph} of the \code{PAFit_net} object. Option \code{plot = "graph"} plots the graph, while \code{plot = "degree"} plots the degree distribution. Option \code{slice} allows selection of the time-step at which the temporal graph is plotted. 
  
  The last two options concern plotting the PA function and node fitnesses (if they are not \code{NULL}).
}
\usage{
\method{plot}{PAFit_net}(x,
     plot = "graph"                         ,
     slice = length(unique(x$graph[,3])) - 1,
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class \code{PAFit_net}.
  }
 \item{plot}{
    String. Possible values are \code{"graph"}, \code{"degree"}, \code{"PA"}, and \code{"fit"}. Default value is \code{"graph"}.
  }
  \item{slice}{
    Integer. Ignored when \code{plot} is not \code{"graph"} or \code{"degree"}. Specifies the time-step at which the graph is plotted. Default value is the final time-step.
  }
    \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}


\value{
  Outputs the desired plot.
}
\author{
  Thong Pham \email{thongpham@thongpham.net}. When \code{plot = "graph"}, the function uses \code{\link{plot.network.default}} in the \pkg{network} package.
}


\examples{
    library("PAFit")
    # a network from Bianconi-Barabasi model
    net        <- generate_BB(N = 50 , m = 10 , s = 10)
    plot(net, plot = "graph")
    plot(net, plot = "degree")
    plot(net, plot = "PA")
    plot(net, plot = "fit")
}
