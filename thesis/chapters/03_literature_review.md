# Literature Review

## Problems we optimized: Crashworthiness models
Crashworthiness models are models of vehicles used to analyze the safety of
their occupants in a crash. @Fig:crashworthiness shows an example. They can be
very elaborate, but very simple ones are also used on initial phases of vehicle
design. These simple ones are one-dimensional Masses-Springs-Dampers Systems
(MSDSs), i.e. they're one-dimensional and are comprised of only ideal masses,
springs and dampers. They can be used to estimate the maximum acceleration a
driver would suffer in a crash. Alternatively, we can look for values of the
springs and dampers that would minimize the maximum acceleration the drive would
suffer in a crash.

![Example of Crashworthiness model. @Tbl:crashworthiness has the legend. Source: [@Mostafa2011-kc]](figs/crashworthiness.png){#fig:crashworthiness width=80%}


|  Mass No. |                    Vehicle components |
|:---------:|--------------------------------------:|
|          1|                    Engine and Radiator|
|          2|             Suspension and Front Rails|
|          3|             Engine Cradle and Shotguns|
|          4| Fire Wall and Part of Body on Its Back|
|          5|                               Occupant|
: Legend of @Fig:crashworthiness. Source: [@Mostafa2011-kc] {#tbl:crashworthiness}

In this work, we set out to analyze the cost/benefits of using MMGA vs NMMGA for
optimization of mechanical systems, and how that changes with respect to the
system's complexity. In order to perform this analysis, we chose to optimize
many systems of different complexities using both algorithms to then analyze the
performance of each algorithm for different system's complexities.
To summarize, the following pseudo-code illustrates our basic workflow:
```{.python caption="Illustration of how we generated the data to compare MMGA and NMMGA"}
for complexity in [1,2,3 ... ]:
  for algorithm in ["MMGA", "NMMGA"]:
    for random_seed in [1,2,3 ...]:
      problem = NewRandomProblem(complexity, random_seed)
      results = Optimize(problem, algorithm)
      SaveResultsToAnalyzeLater(complexity, algorithm, results)
```

Inspired by crashworthiness models, we chose that **the problems we would
optimize were those of finding values of springs and dampers from MSDSs that
minimize the maximum acceleration that a mass - the one farthest away from
impact - would experience with an impact**. For example, for the system of
@Fig:crashworthiness, the optimization problem is: given the values of the
masses and an initial speed of the system, find the values of $k_i$ and $c_i$
that minimize the maximum acceleration $m_5$ will experience with the impact of
the system with the fixed wall on the left.

We chose to optimize only these kind of problems because:

- As explained in @sec:motivation, our focus is not to study the algorithms applied to a specific problem, but to study the algorithms themselves.
- Implementing a solver to this kind of problem is much simpler than, for example, a FEM solver.
- The physics equations which describe these systems, which were used in the PINN metamodels, are relatively simple.
- It's easy to create a problem with arbitrary complexity. @Fig:MsdsExamples shows examples of why that is true: to increase a system's complexity, we can always add more masses/springs/dampers. Also, for a given configuration of springs and dampers, we can easily create infinite amount of problems by picking random mass values.
- By creating a problem with enough complexity, we can have systems which are highly non-linear.

![Examples of MSDSs of increasing complexity (from top to bottom). Source: Author](figs/MsdsExamples.png){#fig:MsdsExamples width=80%}

## Discrete element method solution to MSDS's
[@Mostafa2011-kc]


$$
\begin{bmatrix}
k_{01} & -k_{01}\\
-k_{01} & k_{01}
\end{bmatrix}
\begin{bmatrix}
x_0\\
x_1\\
\end{bmatrix}
+
\begin{bmatrix}
b_{01} & -b_{01}\\
-b_{01} & b_{01}
\end{bmatrix}
\begin{bmatrix}
\dot{x_0}\\
\dot{x_1}\\
\end{bmatrix}
=
\begin{bmatrix}
\ddot{x_0}\\
\ddot{x_1}\\
\end{bmatrix}
$$
{#eq:springLocalK}

@Eq:springLocalK shows test

## Genetic Algorithm

## Meta-models for mechanical optimization
[@Driemeier_undated-za]
[@Gu2018-uk]
[@Gu2018-tf]
[@Lee2022-uz]
[@Wilt2020-np]

## Physics Informed Neural Networks (PINNs)
[@Nascimento2020-xp]


