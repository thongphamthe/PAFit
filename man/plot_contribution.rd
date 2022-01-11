\name{plot_contribution}
\alias{plot_contribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Plotting contributions calculated from the observed data and contributions calculated from simulated data
}
\description{
  This function extracts from a \code{Simulated_Data_From_Fitted_Model} object contributions of rich-get-richer and fit-get-richer effects calculated using simulated networks and plots these contributions versus the contributions calculated from the original observed network. See \code{\link{joint_estimate}} for a description of how the contributions are calculated.
}
\usage{
plot_contribution(simulated_object,
                  original_result,
                  which_plot = "PA",
                  y_label = ifelse("PA" == which_plot,
                  "Contribution of the rich-get-richer effect",
                  "Contribution of the fit-get-richer effect"),
                  legend_pos_x = 0.75,
                  legend_pos_y = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simulated_object}{
    an object of class \code{Simulated_Data_From_Fitted_Model} that contains simulated data.
  }
   \item{original_result}{
    an object of class \code{Full_PAFit_result} that contains the estimation results from the original observed data.
  }
  \item{which_plot}{
    String. ``PA": plots contributions of rich-get-richer effect, ``fit": plots contribution of fit-get-richer effect. Default is ``PA".
  }
   \item{y_label}{
   String. The label for y-axis. Default is "Contribution of rich-get-richer effect". 
  }
\item{legend_pos_x}{Numeric. The horizontal position, between (0,1), of the legend. Default value is \code{0.75}.}  
\item{legend_pos_y}{Numeric. The vertical position, between (0,1), of the legend. Default value is \code{0.9}.} 
}

\value{
 Output a plot.
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
 \code{\link{joint_estimate}}, \code{\link{plot_contribution}}
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

