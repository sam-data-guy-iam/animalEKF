\name{cdlm_robot_twostate_2D}
\alias{cdlm_robot_twostate_2D}
\title{
Shiny app for simulation of 2D robot movement with CDLM and two states.
}
\description{
Shiny app for simulation of 2D robot movement with CDLM and two states.
}
\usage{cdlm_robot_twostate_2D()}
\details{
See \code{cdlm_robot} and \code{cdlm_robot_twostate} for explanation of the basic concepts.  This function simulates a 2-D moving robot with two behavioral states (1 and 2, "slow"/"fast") to model, as well as the switching probabilities between them.

The means of the log-speeds of the two behaviors are simulated by a normal distribution with two means alpha ("unknown true mean of log-speed", types 1 and 2).  The variance in each case is the same and known, as before. The prior means and variances of the velocities are assigned as before. 

The transition probabilities between the behaviors are given by "transition probability between type 1 and 2" and "2 and 1".  If box "are transition probabilities known?" is checked, then they are known.  Otherwise, the transition probabilities will be estimated by a Dirichlet prior (vector "Dirichlet prior values" of form 1->1, 1->2, 2->1, 2->2 of positive numbers, which should roughly correspond to the true probabilities in ratio).  Note that the predictions in this simulation are unlikely to be as good as in the prior 1-D example since there are more parameters to learn and only a limited number of timesteps or particles.

Panel 1 shows the particles' distributions of the velocities for each behavior.  Panel 2 shows the location predictions (black dot and rectangle for type 1, x for type 2 as the widths of the 95\% confidence intervals; again these are of constant witdth since the variance is known). Panel 3 shows the overall resampling weights for the particles.  Panel 4 shows the resampled particles, along with their prediction of location and behavior type. 

Panel 5 shows the behavior-specific particle weights, reflecting howw well the particle's prediction at each behavior matches what was observed.  The overall weights in panel 3 are calculated for each particle as the sum of the behavior-specific weights.  Ideally, the weights in one panel should be high while in the other panel they should be low.  The behavior-specific weights, however, involve the transition probabilities and the previous simulated behavior, so a particle can have a low weight for the correct behavior if it simply predicted the wrong behavior in the previous step.

Panel 6 shows the accuracy of particle predictions of the behavior.  The color (1=black, 2=red) is the true behavior type, and the height of the bar is the fraction of particles correctly predicting it.  Ideally, all bars should be high.

Panel 7 shows convergence over time of the means of the particle distributions of velocity to the true value (vertical dashed line). Ideally these should converge to the true value. 

Panel 8 shows the predicted locations over time (Xs and dots, for each behavior).  Ideally these should gradually line up with the true locations (vertical dashed lines).

If we choose to have the transition probabilities be unknown, Panel 9 shows the distributions (beta is same as one component of a bivariate Dirichlet) of the transition probabilities.  Ideally these should likewise get closer to the true values, the vertical dashed line. 
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Ackerman, Samuel.  "A Probabilistic Characterization of Shark Movement Using Location Tracking Data."  Temple Univeristy doctoral thesis, 2018.  \url{https://digital.library.temple.edu/digital/collection/p245801coll10/id/499150}

Carvalho, Carlos M., Johannes, Michael S., Lopes, Hedibert F., and Nicholas G. Polson.  "Particle learning and smoothing."  Statistical Science, 2010.
}
\note{
Video explanation of simulation applet by author: \url{https://youtu.be/4XR8eB89z7E}

}

\keyword{datasets}
