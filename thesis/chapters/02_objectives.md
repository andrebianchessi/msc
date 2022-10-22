# Objectives
This work sets out to:

1. Define a few mechanical systems with increasing complexity that are comprized
of ideal masses, linear springs and linear dampers.
2. For each one of the mechanical systems, find the parameters (constants of
springs and dampers) that minimize the maximum acceleration that one of the
masses - the one farthest from impact - experiences on an impact with a given
speed. The optimization is done twice for each system using Genetic Algorithms:
first by performing Explicit Time Integration (ETE) of equations found with a
Discrete Element Method (DEM) and, secondly, by using a PINN (Physics Informed
Neural Network) metamodel that approximates the system's dynamic response. The
PINN metamodels are trained without labeled (ground truth) data; the physics
equations are embedded into their training.
3. Redo the optimizations using the same metamodels, without retraining them,
but for different impact speeds (both lower and higher).
4. Analyze and compare the results obtained.

Some of the questions we want to answer are:

- For each combination of *method to find the system's dynamic response (ETE and
  metamodel)*, *system's complexity* and *speed of impact (same, lower than, and
  higher than the one used for metamodel training)*, how did the optimal solution
  found (constants of springs and dampers) and optimal result found (maximum
  acceleration) vary?
- How big was the added cost of training the metamodels? Was it *worth it*, in
  all cases, to use the metamodels; or did the training overhead was greater
  than the efficiency gains for some cases?
- How well did the metamodels perform, i.e. how close to the ETE results were
  their predictions, for the impacts with the speed they were trained on? How
  about for the impacts with different speeds?
