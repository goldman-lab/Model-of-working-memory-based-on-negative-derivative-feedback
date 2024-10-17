# Model-of-working-memory-based-on-negative-derivative-feedback

This is the readme for simulation code associated with the paper:

Lim S, Goldman MS (2013) Balanced cortical microcircuitry for
maintaining information in working memory, Nature
Neuroscience 16:1360-1314

Matlab simulation code written by Sukbin Lim and posted to ModelDB on
8/12/2014

Contents:

FiringRateModel_PM.m: An m-file simulates firing rate models of
parametric working memory circuits based on negative derivative
feedback.  If flag = 1, the external input is transient, resulting a
jump in the network activity and if flag = 0, the external input is
like a step, resulting in ramp-like activity.  The simulations of the
firing rate models were run with a fourth-order explicit Runge-Kutta
method using the function ode45.

SpikingNetworkModel_PM.m: An m-file simulates spiking network models
of parametric working memory circuits based on negative derivative
feedback.  It reproduces Figure 5 in the paper.  The numerical
integration of the network simulation was performed using the
second-order Runge-Kutta algorithm and spike times were approximated
by linear interpolation (Hansel et al. 1998)

- Caution: Simulation of spiking network model in Matlab is slow for a
  large size of network which is required for strong corrective
  feedback.
