# Introduction
TODO:

- Define metamodel (AI that represents a FEM model)
- Show literature examples of metamodels used for genetic algorithm optimization
- Explain that metamodels require a high cost for their training, but might
accelerate the overall optimization for, once trained, the AI is faster than
the FEM model
- We set out to study a mechanical optimization problem using genetic algorithm
(lumped crashworthiness models) but focusing on also analyzing the cost of:
    1. Do the optimization using FEM models
    2. First training a metamodel and then do the simulations using it, instead
    of the FEM model
    3. **Propose** a a hybrid strategy, which starts the problem solution using
    the FEM model, trains the metamodel as data gets generated and, once the
    metamodel is good enough, starts using it for the next genetic algorithm
    generations
- We hypothesize that PINN might be excellent for this case for it can be
trained with a very small dataset
- We also propose to train an AI that replaces the optimization as a whole,
so that the mechanical system at hand can be quickly optimized for other cases