# Motivation {#sec:motivation}

Optimization algorithms are widely spread in all areas of modern Engineering
industry. They allow for reduction of costs and increase of the efficiency and
effectiveness of solutions. The exponential growth in computational power in the
recent years has made it feasible for engineers to optimize problems they could
not optimize before. Still, the complexity of the problems we tackle increases
even faster than the available computational power. Thus, advancements in
optimization techniques are also mandated. In this work, we turn our attention
to optimization by Genetic Algorithms (GAs) applied to mechanical systems. These
algorithms are based on principles of genetics and evolution, and are well
suited for the optimization of problems which deal with a large amount of
variables and have multiple local minima.

However, a caveat GA has is that it needs to constantly evaluate the *fitness*
(i.e. ta score of how good each possible solution is) of each candidate solution 
in the population (i.e. the set of possible solutions) at every iteration. For
example: when using GA to minimize the maximum displacement a structure suffers
at a specific load condition, it's necessary to calculate the structure's
displacement multiple times each iteration of the algorithm. Thus, for expensive
fitness functions, the algorithm's computational cost is proportional to the fitness
function's cost. This might make the algorithm too expensive if the fitness
is computed, for example, using an expensive Finite Element Model (FEM)
simulation.

With this in mind, many recent studies have researched using Machine Learning
(ML) models in conjunction with GA to increase the efficiency of optimizations
of problems which have costly fitness functions. The basic idea is to train a model
which approximates the fitness function but has a much smaller computational cost;
and use it to calculate the fitness in the GA, instead of the original fitness
function. We call these ML models *metamodels* because they are used as
approximations to another model, such as a FEM simulation. There are many very
successful cases in the literature [@Lee2022-uz; @Gu2012-ru; @Gu2018-uk;
@Gu2018-tf; @Driemeier_undated-za].

The challenge of metamodel-based-GAs is that although metamodels can be used for cheaper
evaluation of the fitness function, using them isn't free because there's also a
cost associated with training them; especially when they're trained with
synthetic data generated through numerical simulations, which can be very costly
depending on the problem. If it takes too long to train the metamodels, the
faster evaluation of the fitness function that they provide might not be worth the
added cost that they bring.

Physics-Based Machine Learning is a promising field of research that has been
quickly growing, and has already shown great results when used in conjunction with GAs
[@Zhang2021-vj; @Zhang2022-jp]. Physics-Informed Machine Learning Models (PIMs)
are Machine Learning models which are trained with physics knowledge embedded
into the loss function. This allows them to be trained without labeled data.
Thus, it's possible to train, for example, a metamodel to a physical problem
without needing to run a numerical simulation multiple times beforehand to
obtain training data. This way, PIMs might be good candidates to GAs applied to
problems which have expensive fitness functions. Also, given that they're trained
without labeled data, PIMs might have good extrapolation capabilities. This
characteristic might make then perform well even on a very large domain of
possible solutions, which is exactly the cases in which GAs excel.

In this work, we set out to analyze how much performance improvement can be
obtained for a GA applied to mechanical systems when PIMs are used;.

Note that our end goal is not to assess the applicability of this technique to a
specific problem, but to further investigate the technique itself to understand
its potentials and limitations. For this reason, we only analyze
Mass-Spring-Damper Systems. Solvers for these systems have a small
implementation overhead, which allows us to spend less effort on implementing
FEM solvers for example. Besides, it's easy to generate arbitrarily
simple/complex systems of this kind.
