# Introduction
Optimization algorithms are widely spread in all areas of modern Engineering industry.
They allow for reduction of costs and increase of efficiency and effectiveness of
solutions. The exponential growth in computational power in the recent years has
made it feasible for engineers to optimize problems they could not optimize
before. Still, the complexity of the problems we tackle increases even
faster than the available computational power. Thus, advancements in optimization
techniques are also mandated.
In this work, we turn our attention to optimization by Genetic Algorithms (GAs)
applied to mechanical systems. *The genetic algorithm is an optimization and search technique based on the principles of genetics and natural selection* [@Haupt2004-hj]. Some of its advantages, when compared to gradient base methods, are [@Haupt2004-hj]:
- Doesn’t require derivative information
- Simultaneously searches from a wide sampling of the cost surface
- Deals with a large number of variables
- Is well suited for parallel computers
- Optimizes variables with extremely complex cost surfaces (they can jump out of local minima)

However, a caveat GA has is that it needs to constantly evaluate the *cost
function* (i.e. the function that returns quality of a possible solution) of the
whole population (i.e. a set of possible solutions) at every iteration. For
example, when using GA to minimize the maximum displacement a structure suffers
at a specific condition, it's necessary to calculate the structure's
displacement multiple times, with Finite Element Method (FEM), for example, at each
iteration of the algorithm. Thus, for expensive cost functions, the algorithm's
computational cost is proportional to the *cost function*'s cost. This might
make the algorithm too expensive if the cost function is found, for example,
using a very refined FEM simulation.

With this in mind, many recent studies have researched using Machine Learning
(ML) models in conjunction with GA to increase the efficiency of optimizations
of problems which have costly cost functions. The basic idea is to train a model
which approximate the cost function but has a much smaller computational cost;
and use it to calculate the cost in the GA, instead of the original cost
function. For cases in which the cost function itself is a model (i.e. and
approximation of the reality) such as a FEMs, the ML models
used to approximate them can be called *metamodels*. There are many very
successful cases in the literature [@Wilt2020-np, @Lee2022-uz, @Gu2012-ru,
@Gu2018-uk, @Gu2018-tf, @Driemeier_undated-za].

What we have not found in the literature, however, are in depth analyses of how
the efficiency of metamodel-based-GA optimization changes with respect to the problem's complexity and the ML model used.
After all, there's also a cost associated with training the metamodel. Can in be that, in some cases, we'd need to calculate the cost function $1$M times to train a metamodel, but would only calculate it $10$k times in the GA? In those cases, using a metamodel would not make sense.

Physics-Based Machine Learning is a promising field of research that has been quickly growing,
and has already shown great results when used in metamodel-based-GA [@Zhang2021-vj, @Zhang2022-jp]. Physics-Informed-Neural-Networks (PINNs) are
Neural Networks (NNs) which have physics knowledge embedded into their loss-function. This allows them to be trained without labeled data. Thus, it's possible to train, for example, a metamodel to an expensive FEM simulation without needing to run the FEM simulation multiple times beforehand. This way, PINNs might be ideal to metamodel-based-GA applied to problems which have expensive cost functions.

Also, given that they're trained without labeled data, PINNs might have better extrapolation capabilities than NNs. This characteristic might make then perform well even on a very large domain of possible solutions, which is exactly the cases in which GAs excel.

In this work, we set out to analyze how much performance improvement can be obtained for a GA applied to mechanical systems when PINN metamodels are used; and how that varies with respect to the system's and metamodel's complexity. We also want to investigate the extrapolation capabilities of PINNs in this technique, i.e. how they perform when evaluating conditions much different than the ones they were trained on.

Note that our focus is not to investigate this technique applied to a specific problem, but of the technique itself. We know that results to some case studies are not trivially generalizable, but they may still increase our understanding on potentials and limitations of this promising new technique, provide some insights and new ideas worth pursuing, and facilitate future work that might be similar to this one.

For this reason, we only analyze Mass-Spring-Damper Systems (MSDSs). Solvers for these systems require a small implementation overhead, which allows us to spend less effort on implementing FEM solvers for example. Besides, it's easy to generate arbitrarily simple/complex systems of this kind.
