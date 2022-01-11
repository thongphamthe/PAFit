\name{generate_simulated_data_from_estimated_model}
\alias{generate_simulated_data_from_estimated_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Generating simulated data from a fitted model  
}
\description{
  This function generates simulated networks from a fitted model and performs estimations on these simulated networks with the same setting used in the original estimation. Each simulated network is generated using parameters of the fitted model, while keeping other aspects of the growth process as faithfully as possible to the original observed network.
}
\usage{
generate_simulated_data_from_estimated_model(net_object, net_stat, result, M = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net_object}{
    an object of class \code{PAFit_net} that contains the original network.
  }
  \item{net_stat}{
    An object of class \code{PAFit_data} which contains summarized statistics of the original network. This object is created by the function \code{\link{get_statistics}}. 
  }
   \item{result}{
    An object of class \code{Full_PAFit_result} which contains the fitted model obtained by applying the function \code{\link{joint_estimate}}. 
  }
\item{M}{integer. The number of simulated networks. Default value is \code{5}.}  
}

\value{
  Outputs a \code{Simulated_Data_From_Fitted_Model} object, which is a list containing the following fields:
  \itemize{
    \item \code{graph_list}: a list containing \code{M} simulated graphs.
    
    \item \code{stats_list}: a list containing \code{M} objects of class \code{PAFit_data}, which are the results of applying \code{\link{get_statistics}} on the simulated graphs.
    \item \code{result_list}: a list containing \code{M} objects of class \code{Full_PAFit_result}, which are the results of applying \code{\link{joint_estimate}} on the simulated graphs.
}
}

\author{
  Thong Pham \email{thongphamthe@gmail.com}
}
\references{
  1. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. (\doi{10.1371/journal.pone.0137796}).
  
  2. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. (\doi{10.1038/srep32558}).
  
 3. Pham, T., Sheridan, P. & Shimodaira, H. (2020). PAFit: An R Package for the Non-Parametric Estimation of Preferential Attachment and Node Fitness in Temporal Complex Networks. Journal of Statistical Software 92 (3). (\doi{10.18637/jss.v092.i03}).

4. Inoue, M., Pham, T. & Shimodaira, H. (2020). Joint Estimation of Non-parametric Transitivity and Preferential Attachment Functions in Scientific Co-authorship Networks. Journal of Informetrics 14(3). (\doi{10.1016/j.joi.2020.101042}).
}
\seealso{
 \code{\link{get_statistics}}, \code{\link{joint_estimate}}, \code{\link{plot_contribution}}
}

\examples{
\dontrun{
  
  library("PAFit")
  net_object     <- generate_net(N = 500, m = 10, s = 10, alpha = 0.5)
  net_stat       <- get_statistics(net_object) 
  result         <- joint_estimate(net_object, net_stat)
  simulated_data <- generate_simulated_data_from_estimated_model(net_object, net_stat, result)
  plot_contribution(simulated_data, result, which_plot = "PA")
  plot_contribution(simulated_data, result, which_plot = "fit")
  }
}

\concept{preferential attachment}
\concept{attachment function}
\concept{fitness}

