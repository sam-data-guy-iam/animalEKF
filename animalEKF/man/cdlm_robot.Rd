\name{cdlm_robot}
\alias{cdlm_robot}
\title{
Shiny app for 1D simulation of robot movement with CDLM.
}
\description{
Shiny app for 1D simulation of robot movement with CDLM.
}
\usage{cdlm_robot()}
\details{
This shiny app illustrates a 1-D robot movement model.  Here, T ("maximum number of iterations") steps are simulated for a robot moving along a 1-D line.  Each of the T steps represents a length of time represented by "time step (sec)" seconds; the longer the interval, the more location uncertainty there is between steps.  

At each step, the robot moves with velocity (v_t) modeled by a normal distribution with mean alpha ("unknown true mean of velocity") and variance "known true variance." For simplicity, we will only attempt to model the mean velocity while sequentially observing only the locations, since we assume the variane is known. 

The particle filter learns the movement parameters through N ("number of particles") particles, or independent simulations.  At each point in time, the filter simulates N draws of the velocity from the prior distribution, a normal distribution with mean mu ("prior mean on velocity mean") and variance sigma ("prior variance on velocity mean"). These distributions are shown in color panel 1.  Ideally, over time the colored distributions should converge to the true one (thick black curve).

Panel 2 shows each particle's prediction of the location (black dot) and the 95\% confidence interval of this prediction (width of colored rectangle--in this case all are the same width because the variance is known).  The true observed location and the previous one are the two dashed vertical lines.  Particles whose dots fall closer to the vertical line have better prediction.

Note: this simulation works best if the distribution of true velocity is either clearly negative or positive.  If the distribution straddles v_t=0 with significant probability, movement will be more difficult to visualize.  In the above case, the location should be the one more in the direction of the sign of velocity (i.e. if velocity distribution >0, then the robot should be consistently moving to the right).

Panel 3 shows the weights of the particles as calculated by the closeness of their location prediction to the true one.  Closer particles in panel 2 should have higher weights (the colors correspond).  In panel 4 the particles are resampled by their weighs (with replacement) and the predictions are shown.  More of the predictions should be closer to the truth here than in panel 2.

Panel 5 shows convergence over time of the means of the particle distributions of velocity to the true value.  Panel 6 shows the guesses over time of the particles as to the true location.  Similarly, we should see that over time these predictions line up with the true locations (vertical dashed lines).

The particle filter models the true value of location and the true velocity.  It is reasonable to assume that these are independent, which is why the covariance matrices given are diagonal. 
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Ackerman, Samuel.  "A Probabilistic Characterization of Shark Movement Using Location Tracking Data."  Temple Univeristy doctoral thesis, 2018.  \url{https://digital.library.temple.edu/digital/collection/p245801coll10/id/499150}

Carvalho, Carlos M., Johannes, Michael S., Lopes, Hedibert F., and Nicholas G. Polson.  "Particle learning and smoothing."  Statistical Science, 2010.

}
\note{
Video explanation of simulation applet by author: \url{https://youtu.be/iVG_bCU0jCA}

}
\examples{
}
\keyword{datasets}
