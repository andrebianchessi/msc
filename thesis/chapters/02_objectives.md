# Objectives
This work sets out to:

1. Define a few mechanical systems with increasing complexity that are comprized
of ideal masses, linear springs and linear dampers.
2. For each one of the mechanical systems, find the parameters (constants of
springs and dampers) that minimize the maximum acceleration that one of the
masses experiences on an impact with a given speed. The
optimization is done twice for each system using a Genetic Algorithm: first
using a DES (discrete element simulation) and, secondly, a PINN (Physics
Informed Neural Network) metamodel to find the system's dynamic response. The
PINN metamodels are trained without labeled data.
3. Redo the optimizations using the same PINNs, without retraining them, but for
different impact speeds (both lower and higher).
4. Analyze and compare the results obtained.

Some of the questions we want to answer are:

- For each combination of *method to find the system's dynamic response (DES and
  PINN)*, *system's complexity* and *speed of impact (same, lower than, and
  higher than the one used for PINN training)*, how did the optimal solution
  found (constants of springs and dampers) and optimal result found (maximum
  acceleration) vary?
- How big was the added cost of training the metamodels? Was it *worth it*, in
  all cases, to use the metamodels; or did the training overhead was greater
  than the efficiency gains for some cases?
- How well did the metamodels perform, i.e. how close to the DES results were
  their predictions, for the impacts with the speed they were trained on? How
  about for the impacts with different speeds?
