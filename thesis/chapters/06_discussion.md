# Discussion {#sec:discussion}

In this work, we explored the potential of Physics-Informed Machine Learning Models (PIMs) to bolster the efficiency of Genetic Algorithms (GAs)
for the optimization of mechanical systems. Traditionally, these GAs use numerical simulations
(such as Explicit Time Integration or the Finite Element Method) to measure the fitness of each candidate solution. These numerical
simulations can be quite expensive and become computational bottlenecks, but - once trained - PIMs can rapidly estimate the fitness
of candidates. The training of these models is, however, not free; so at first glance it's not easy to say if using them
is indeed efficient or not.

Our investigations revealed many cases in which *PIM-enhanced GA* yielded substantial time savings when compared to the
*Numerical GA* (which uses Explicit Time Integration to evaluate each candidate solution).
Moreover, we saw that *PIM-enhanced GAs* were not only faster but also capable of producing solutions of comparable quality, or even superior quality, to those obtained using the *Numerical GA*.

In some of the experiments, however, *PIM-enhanced GA* led to a solution which was worse than a random guess;
so it's not always that this method works. Also, finding the hyperparameters which cause the models to train fast enough while
still causing them to be good-enough estimators is tricky and very operator-dependent. Thus, we can't argue that this is
a technique that should always replace *Numerical GA*.

However, it's vital to recognize that these findings do not diminish the promise of *PIM-enhanced GAs*. After all, 11 out of 15 experiments
showed not only a much smaller execution time but
a larger *Mean Score* (the mean between the efficiency and the quality score) for the *PIM-enhanced GA*, and in one
case the solution found by the *PIM-enhanced GA* was even better than the one found by the *Numerical GA*.
Some clear follow-up questions worthy of investigation are:

- *How does the performance of the technique changes with more refined models?* We only used simple linear models, so it stands to question if more refined models can be trained faster but reach the same inference quality.

- *How does it perform in other mechanical problems?* We only analyzed systems made up of masses springs and dampers, but the technique can easily be applied to other problems as well.

- *What are the optimal hyperparameters for the training? How do we find them for each specific problem?*

In conclusion, our work establishes a compelling case for further research about the incorporation of
Physics-Informed Machine Learning Models into GAs to dramatically enhance optimization performance. This approach not only holds the potential to significantly accelerate solution searches but also to maintain or even improve solution quality, which can unlock possibilities for tackling previously intractable problems across a multitude of disciplines.