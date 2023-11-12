# Objectives {#sec:objectives}
This work sets out to:

1. Define a few mechanical systems with different complexities that are
comprized of ideal masses, linear springs and linear dampers.
2. For each one of the mechanical systems, find the parameters (constants of
springs and dampers) that minimize the maximum acceleration that one of the
masses - the one farthest from impact - experiences on an impact with a given
speed. The optimization is done twice for each system using Genetic Algorithms:
first by performing Explicit Time Integration (ETE) of equations found with a
Discrete Element Method (DEM) and, secondly, by using a PIM (Physics-Informed
Machine Learning Model) metamodel that approximates the system's dynamic
response. The first optimization method is referred to as *P-GA*, and the second
one is referred to as *E-GA*. 
3. Analyze and compare the results obtained.

Some of the questions we want to answer are:

- How big was the added cost of training the metamodels? Was it *worth it*, in
  all cases, to use the metamodels; or did the training overhead was greater
  than the efficiency gains for some cases?
- How well did the metamodels perform, i.e. how close to the ETE results were
  their predictions, for the impacts with the speed they were trained on? How
  about for the impacts with different speeds?
- How does the system complexity impact the efficiency of using PIMs? I.e. do
  they provide a bigger efficiency gain when the system complexity is larger?
