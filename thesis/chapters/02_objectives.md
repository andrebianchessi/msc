# Objectives {#sec:objectives}
This work sets out to:

1. Perform case studies of *dynamic mechanical system optimization* using
two approaches: *P-GA* and *E-GA*. Both are genetic algorithms, but while *E-GA*
evaluates the fitness of the solutions using explicit time integration ([ETI](#sec:eti)),
*P-GA* does so by using PIMs (Physics Informed Machine Learning Models) that
describe the time response of the system as a function of *time* and of the
system's properties. 
2. Compare the performance and the results obtained by each method.

Some of the questions we seek to answer are:

- How well do the PIMs perform when compared to ETI (Explicit time integration)?
Are they good approximators of the system's time response?
- How big was the added cost of training the models for P-GA? Is the added cost
of training the models worth the faster evaluation time that they provide?
- How good are the solutions found with P-GA when compared to E-GA?
