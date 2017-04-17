\name{NEWS}
\alias{NEWS}
\title{Package News}

\section{version 0.9.8.8 (Release date: 2017-04-17)}{
  \subsection{NEW}{
    \itemize{
    \item a new interface for almost every main functionality of the package
    \item a completely re-written tutorial 
    }
  }
}


\section{version 0.9.7.5 (Release date: 2017-02-15)}{
  \subsection{FIXED}{
    \itemize{
    \item{a bug in plotting the PA function when the degree threshold is set too high} 
    }
  }
    \subsection{CHANGED}{
    \itemize{
      \item{improved the exposition in the documentation}
      \item{improved code readability}
    }
  }
  \subsection{NEW}{
    \itemize{
    \item wrapper functions to generate networks from well-known models 
    \item complete overhaul of the main functions. 
    }
  }
}


\section{version 0.9.6 (Release date: 2017-02-02)}{
  \subsection{FIXED}{
    \itemize{
      \item{minor bugs in edge cases}
    }
  }
    \subsection{CHANGED}{
    \itemize{
      \item{improved the exposition in the documentation}
      \item{improved code readability}
    }
  }
  \subsection{NEW}{
    \itemize{
      \item{adding confidence interval for the estimated attachment exponenent in all methods}
    }
  }
}


\section{version 0.9.5 (Release date: 2017-02-01)}{
  \subsection{FIXED}{
    \itemize{
      \item{some minor bugs}
    }
  }
    \subsection{CHANGED}{
    \itemize{
      \item{improved the exposition in the documentation}
      \item{improved code readability}
    }
  }
  \subsection{NEW}{
    \itemize{
      \item{starting to use GitHub for bugs reporting: \url{https://github.com/thongphamthe/PAFit/issues}}
      \item{starting to use a parallel GitHub homepage: \url{https://github.com/thongphamthe/PAFit}}
      \item{the corrected Newman's method and Jeong's method are implemented}
      \item{functions to evaluate error between estimation and true values}
    }
  }
}


\section{version 0.9.3 (Release date: 2017-1-20)}{
  \subsection{FIXED}{
    \itemize{
      \item{a bug in reading data in which time-stamps are factors}
    }
  }
  
  \subsection{NEW}{
    \itemize{
      \item{functions for performing cross validation when assuming \eqn{A_k = k^\alpha}. This will reduce the cross validation to one-dimensional search, since in this case we only have to find the optimum parameter \eqn{s} for fitness}
    \item{starting use this NEWS}
  }
}
}

