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

### Base case: Local Matrices
@Fig:discreteElementSimple shows the simplest system that contains ideal masses,
springs and dampers. The spring has its natural/relaxed length when $x_0=x_1=0$.

![MSDS of 2 masses, 1 spring and 1 damper. Source: Author](figs/discreteElementSimple.png){#fig:discreteElementSimple width=80%}

From Newton's second law, we obtain @eq:demLocalEq, in which $\dot{x_i}$ and $\ddot{x_i}$
represent $x_i$'s first and second time derivative - speed and acceleration - respectively.

$$
\begin{bmatrix}
-k_{01} & k_{01}\\
k_{01} & -k_{01}
\end{bmatrix}
\begin{bmatrix}
x_0\\
x_1\\
\end{bmatrix}
+
\begin{bmatrix}
-c_{01} & c_{01}\\
c_{01} & -c_{01}
\end{bmatrix}
\begin{bmatrix}
\dot{x_0}\\
\dot{x_1}\\
\end{bmatrix}
=
\begin{bmatrix}
m_0 & 0\\
0 & m_1
\end{bmatrix}
\begin{bmatrix}
\ddot{x_0}\\
\ddot{x_1}\\
\end{bmatrix}
$$
{#eq:demLocalEq}

### Induction: Assembling global matrices

By adding another mass, spring and damper, we can obtain the system at @fig:discreteElementSimple2.

![MSDS of 2 masses, 2 spring and 2 dampers. Source: Author](figs/discreteElementSimple2.png){#fig:discreteElementSimple2 width=80%}

Equations of motion of system illustrated at @fig:discreteElementSimple2, also obtained by Neuton's second law, are:

$$
\begin{aligned}
&
\begin{bmatrix}
-k_{01} &  k_{01} & 0\\
k_{01}  & -k_{01} {\color{blue}- k_{12}} & {\color{blue} k_{12}} \\
0       &  {\color{blue}k_{12}}          & {\color{blue}-k_{12}} \\
\end{bmatrix}
\begin{bmatrix}
x_0\\
{\color{blue}x_1}\\
{\color{blue}x_2}
\end{bmatrix} \\
& +
\begin{bmatrix}
-c_{01} &  c_{01} & 0\\
c_{01}  & -c_{01} {\color{blue}- c_{12}} & {\color{blue} c_{12}} \\
0       &  {\color{blue}c_{12}}          & {\color{blue}-c_{12}} \\
\end{bmatrix}
\begin{bmatrix}
\dot{x_0}\\
{\color{blue}\dot{x_1}}\\
{\color{blue}\dot{x_2}}
\end{bmatrix}
= \\
&
\begin{bmatrix}
m_0 & 0 &   0\\
0 & {\color{blue}m_1} &   0\\
0 & 0   & {\color{blue}m_2}\\
\end{bmatrix}
\begin{bmatrix}
\ddot{x_0}\\
{\color{blue}\ddot{x_1}}\\
{\color{blue}\ddot{x_2}}
\end{bmatrix}
\end{aligned}
$$
{#eq:dem2Masses}

By adding another mass, spring and damper connected to $m_1$, we get the system at @fig:discreteElementSimple3.

![MSDS of 3 masses, 3 spring and 3 dampers. Source: Author](figs/discreteElementSimple3.png){#fig:discreteElementSimple3 width=80%}

The equations of motion of the system illustrated at @fig:discreteElementSimple3 are:

$$
\begin{aligned}
\begin{bmatrix}
-k_{01} &  k_{01}                    &        0&          0\\
k_{01}  & -k_{01} - k_{12} {\color{blue}- k_{13}}  & k_{12}  &   {\color{blue}  k_{13}}\\
0       &  k_{12}                    & -k_{12} &          0\\
0       &                    {\color{blue}  k_{13}}&        0&    {\color{blue} -k_{13}}\\
\end{bmatrix}
\begin{bmatrix}
x_0\\
{\color{blue}x_1}\\
x_2\\
{\color{blue}x_3}
\end{bmatrix}
+ & \\
\begin{bmatrix}
-c_{01} &  c_{01}                    &        0&          0\\
c_{01}  & -c_{01} - c_{12} {\color{blue}- c_{13}}  & c_{12}  &   {\color{blue}  c_{13}}\\
0       &  c_{12}                    & -c_{12} &          0\\
0       &                    {\color{blue}  c_{13}}&        0&    {\color{blue} -c_{13}}\\
\end{bmatrix}
\begin{bmatrix}
\dot{x_0}\\
{\color{blue}\dot{x_1}}\\
\dot{x_2}\\
{\color{blue}\dot{x_3}}
\end{bmatrix}
=&\\
\begin{bmatrix}
m_0 & 0 &   0 & 0  \\
0 & {\color{blue}m_1} &   0 & 0  \\
0 & 0   & m_2 & 0  \\
0 & 0   &   0 & {\color{blue}m_3}
\end{bmatrix}
&
\begin{bmatrix}
\ddot{x_0}\\
{\color{blue}\ddot{x_1}}\\
\ddot{x_2}\\
{\color{blue}\ddot{x_3}}
\end{bmatrix}
\end{aligned}
$$
{#eq:dem3Masses}

By adding just an extra spring and a damper connecting $m_0$ and $m_3$, we get the system at @fig:discreteElementSimple4.

![MSDS of 3 masses, 4 spring and 4 dampers. Source: Author](figs/discreteElementSimple4.png){#fig:discreteElementSimple4 width=80%}

The equations of motion of the system illustrated at @fig:discreteElementSimple4 are:

$$
\begin{aligned}
\begin{bmatrix}
-k_{01} {\color{blue}-k_{03}} &  k_{01}                    &        0&             {\color{blue}k_{03}}\\
k_{01}          & -k_{01} - k_{12} - k_{13}  & k_{12}  &             k_{13}\\
0               &  k_{12}                    & -k_{12} &                  0\\
{\color{blue}k_{03}}          &                      k_{13}&        0&     -k_{13}{\color{blue}-k_{03}}\\
\end{bmatrix}
\begin{bmatrix}
{\color{blue}x_0}\\
x_1\\
x_2\\
{\color{blue}x_3}
\end{bmatrix}
+ & \\
\begin{bmatrix}
-c_{01} {\color{blue}-c_{03}} &  c_{01}                    &        0&             {\color{blue}c_{03}}\\
c_{01}          & -c_{01} - c_{12} - c_{13}  & c_{12}  &             c_{13}\\
0               &  c_{12}                    & -c_{12} &                  0\\
{\color{blue}c_{03}}          &                      c_{13}&        0&     -c_{13}{\color{blue}-c_{03}}\\
\end{bmatrix}
\begin{bmatrix}
{\color{blue}\dot{x_0}}\\
\dot{x_1}\\
\dot{x_2}\\
{\color{blue}\dot{x_3}}
\end{bmatrix}
=&\\
\begin{bmatrix}
{\color{blue}m_0} & 0 &   0 & 0  \\
0 & m_1 &   0 & 0  \\
0 & 0   & m_2 & 0  \\
0 & 0   &   0 & {\color{blue}m_3}
\end{bmatrix}
&
\begin{bmatrix}
{\color{blue}\ddot{x_0}}\\
\ddot{x_1}\\
\ddot{x_2}\\
{\color{blue}\ddot{x_3}}
\end{bmatrix}
\end{aligned}
$$
{#eq:dem2Masses2}

## Genetic Algorithm

## Meta-models for mechanical optimization
[@Driemeier_undated-za]
[@Gu2018-uk]
[@Gu2018-tf]
[@Lee2022-uz]
[@Wilt2020-np]

## Physics Informed Neural Networks (PINNs)
[@Nascimento2020-xp]


