\name{PAFit_oneshot}
\alias{PAFit_oneshot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Estimating the nonparametric preferential attachment function from one single snapshot.  
}
\description{
  This function estimates the attachment function \eqn{A_k} from one snapshot.
}
\usage{

PAFit_oneshot(net  , 
              M    = 100,
              S    = 5,
              loop = 5,
              G = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
    a net.
  }
   \item{M}{
  Integer. Number of simulated networks in each iteration. Default is \code{100}.
  }
  \item{S}{
  Integer. Number of iterations inside each loop. Default is \code{5}.
  }
  \item{loop}{
  Integer. Number of loops of the whole process. Default is \code{5}.
  }
    \item{G}{
  Integer. Number of bins for the PA function. Default is \code{1000}.
  }
}

\value{
  Outputs the estimated result.
}

\author{
  Thong Pham \email{thongphamthe@gmail.com}
}
%\references{
%  1. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).
  
%  2. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{http://www.nature.com/articles/srep32558}).
  
% 3. Pham, T., Sheridan, P. & Shimodaira, H. (2020). PAFit: An R Package for the Non-Parametric Estimation of Preferential Attachment and Node Fitness in Temporal Complex Networks. Journal of Statistical Software 92 (3), doi:10.18637/jss.v092.i03. (\url{http://dx.doi.org/10.18637/jss.v092.i03})
%}
%\seealso{
%  See \code{\link{get_statistics}} for how to create summerized statistics needed in this function.
  
%  See \code{\link{Jeong}}, \code{\link{Newman}} and \code{\link{only_A_estimate}} for functions to estimate the attachment function in isolation.
  
%    See \code{\link{only_F_estimate}} for a function to estimate node fitnesses in isolation.
%}

\examples{
\dontrun{
  library("PAFit")
  }
}

\concept{preferential attachment}
\concept{attachment function}

