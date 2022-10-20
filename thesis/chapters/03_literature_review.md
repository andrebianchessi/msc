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

### Local Matrices
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

Rewriting the above equation:
$$
K_{l}
\begin{bmatrix}
x_0\\
x_1\\
\end{bmatrix}
+
C_{l}
\begin{bmatrix}
\dot{x_0}\\
\dot{x_1}\\
\end{bmatrix}
=
M_{l}
\begin{bmatrix}
\ddot{x_0}\\
\ddot{x_1}\\
\end{bmatrix}
$$
{#eq:demLocalMatrices}

The matrices $K_{l}$, $C_{l}$ and $M_{l}$ are termed, respectively, **local stiffness matrix**, **local damping matrix** and **local mass matrix**.

### Assembling global matrices

By adding another mass, spring and damper to @fig:discreteElementSimple, we can obtain the system at @fig:discreteElementSimple2.

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
{#eq:dem1}

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
{#eq:dem2}

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
{#eq:dem3}

By looking at [@eq:dem1; @eq:dem2; @eq:dem3], we can see that the equations of
motion of the whole system are obtained by superposing the local matrices of
each element, shown at @eq:demLocalEq. The local matrices that are being added
are highlighted in $\textcolor{blue}{\text{blue}}$.

This is the same process done at the Finite Element Method (FEM) to assemble the
global matrices from the elements' local matrices. The reader can read in more
detail at [@Logan2007-bq; @Alves2020-gz].

### Dynamic simulation: ODE
One way to finding the dynamic response of a system is to find the system's Ordinary Differential Equation (ODE), i.e. a system's state vector $X$ and an expression to calculate the state vector's time derivative $\dot{X}$. If we have initial values of the state vector and know how to calculate the its time derivative, we can use numerical methods such as Forward Euler or Runge-Kutta to integrate the ODE, i.e. find approximate values of $X$ after the initial time instant.

We can represent the state of the system with the state vector $X$:
$$
X =
\begin{bmatrix}
x_0\\
x_1\\
\vdots\\
x_n\\
\dot{x_0}\\
\dot{x_1}\\
\vdots\\
\dot{x_n}
\end{bmatrix}
$$

It's time derivative is given by:
$$
\dot{X} =
\begin{bmatrix}
\dot{x_0}\\
\dot{x_1}\\
\vdots\\
\dot{x_n}\\
\ddot{x_0}\\
\ddot{x_1}\\
\vdots\\
\ddot{x_n}
\end{bmatrix}
$$

Notice that the first half of $\dot{X}$ is, simply, the second half of the $X$. Hence, we just need to find expressions for the second derivatives of the displacements ($\ddot{x_1}$, ... , $\ddot{x_n}$).

After assembling the global matrices, we'll be left with a matrix equation in the
following form:
$$
K
\begin{bmatrix}
x_0\\
x_1\\
\vdots\\
x_n
\end{bmatrix}
+
C
\begin{bmatrix}
\dot{x_0}\\
\dot{x_1}\\
\vdots\\
\dot{x_n}
\end{bmatrix}
=
M
\begin{bmatrix}
\ddot{x_0}\\
\ddot{x_1}\\
\vdots\\
\ddot{x_n}
\end{bmatrix}
$$

We can always multiply both sides by $M^{-1}$, since $M$ is a diagonal matrix with strictly positive numbers in the diagonal. This way, the system's ODE is defined by the two following equations:
$$
X =
\begin{bmatrix}
x_0\\
x_1\\
\vdots\\
x_n\\
\dot{x_0}\\
\dot{x_1}\\
\vdots\\
\dot{x_n}
\end{bmatrix}
\quad\quad\text{,}\quad\quad
\dot{X} =
\begin{bmatrix}
\dot{x_0}\\
\dot{x_1}\\
\vdots\\
\dot{x_n}\\
\ddot{x_0}\\
\ddot{x_1}\\
\vdots\\
\ddot{x_n}
\end{bmatrix}
\quad\quad\text{ and }\quad\quad
\begin{bmatrix}
\ddot{x_0}\\
\ddot{x_1}\\
\vdots\\
\ddot{x_n}
\end{bmatrix}
=
M^{-1}
\left(
K
\begin{bmatrix}
x_0\\
x_1\\
\vdots\\
x_n
\end{bmatrix}
+
C
\begin{bmatrix}
\dot{x_0}\\
\dot{x_1}\\
\vdots\\
\dot{x_n}
\end{bmatrix}
\right)
$$
{#eq:ode}

Note that inverting $M$ is trivial, since it's diagonal.

The boundary conditions can then be applied directly at the $\dot{X}$ vector. For fixed masses, for example, we just need to replace the appropriate element of $\dot{X}$ with zero.

#### Example {.unnumbered}
Taking the system of @fig:discreteElementSimple. Considering that the mass on
the left is fixed and the one at the right starts with a positive initial
displacement:
$$
\begin{aligned}
&m_0 \text{ is fixed} \\
&k_{01} = c_{01} = 1 \\
&m_1=1 \\
&x_0|_{t=0} = x_1|_{t=0} = 10 \\
&\dot{x_1}|_{t=0} = 0 
\end{aligned}
$$
{#eq:edoInit}

The state vector for this problem is:
$$
X =
\begin{bmatrix}
x_0\\
x_1\\
\dot{x_0}\\
\dot{x_1}
\end{bmatrix}
$$

By assembling the system's global matrices (which was already done at @eq:demLocalEq), from @eq:ode we have:
$$
\begin{bmatrix}
\ddot{x_0}\\
\ddot{x_1}\\
\end{bmatrix}
=
\begin{bmatrix}
1/m_0 & 0\\
0 & 1/m_1
\end{bmatrix}
\left(
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
\right)
$$

Replacing the values from @eq:edoInit:
$$
\begin{bmatrix}
\ddot{x_0}\\
\ddot{x_1}\\
\end{bmatrix}
=
\begin{bmatrix}
1/m_0& 0\\
0 & 1
\end{bmatrix}
\left(
\begin{bmatrix}
-1 & 1\\
1 & -1
\end{bmatrix}
\begin{bmatrix}
x_0\\
x_1\\
\end{bmatrix}
+
\begin{bmatrix}
-1 & 1\\
1 & -1
\end{bmatrix}
\begin{bmatrix}
\dot{x_0}\\
\dot{x_1}\\
\end{bmatrix}
\right)
$$

$$
\begin{bmatrix}
\ddot{x_0}\\
\ddot{x_1}\\
\end{bmatrix}
=
\begin{bmatrix}
1/m_0(-x_0+x_1-\dot{x_0}+\dot{x_1})\\
x_0-x_1+\dot{x_0}-\dot{x_1}
\end{bmatrix}
$$

Since we considered $m_0$ to be fixed, we replace the expression of $\ddot{x_0}$
with $0$:

$$
\begin{bmatrix}
\ddot{x_0}\\
\ddot{x_1}\\
\end{bmatrix}
=
\begin{bmatrix}
0\\
x_0-x_1+\dot{x_0}-\dot{x_1}
\end{bmatrix}
$$

Thus, the EDO is defined as:
$$
X =
\begin{bmatrix}
x_0\\
x_1\\
\dot{x_0}\\
\dot{x_1}
\end{bmatrix}
\quad\quad\text{ and }\quad\quad
\dot{X} =
\begin{bmatrix}
\dot{x_0}\\
\dot{x_1}\\
0\\
x_0-x_1+\dot{x_0}-\dot{x_1}
\end{bmatrix}
$$

If the EDO were to be integrated with a simple Forward Euler using a $0.1$ time-step, the first time-step would be:

- **t = 0** (initial values taken from @eq:edoInit)

$$
X|_{t=0} =
\begin{bmatrix}
0\\
10\\
0\\
0
\end{bmatrix}
$$

$$
\dot{X} =
\begin{bmatrix}
\dot{x_0}\\
\dot{x_1}\\
0\\
x_0-x_1+\dot{x_0}-\dot{x_1}
\end{bmatrix}
\rightarrow
\dot{X}|_{t_0} =
\begin{bmatrix}
0\\
0\\
0\\
-10
\end{bmatrix}
$$

- **t = 0.1**

$$
X|_{t=0.1} =
X|_{t=0} + 0.1 \cdot \dot{X}|_{t_0} = 
\begin{bmatrix}
0\\
10\\
0\\
0
\end{bmatrix}
+ 0.1
\begin{bmatrix}
0\\
0\\
0\\
-10
\end{bmatrix}
$$

$$
X|_{t=0.1} =
\begin{bmatrix}
0\\
10\\
0\\
-1
\end{bmatrix}
$$

### Usage in code

At the code's [repository](https://github.com/andrebianchessi/msc), the file
```software/problem_test.cc``` contains many examples of how the software we
implemented can be used. @Mostafa2011-kc was extensively used as reference for
implementing test cases.

Basically, we first initialize a `Problem` object. Then, we use the `AddMass`,
`AddSpring` and `AddDampers` methods to add masses, springs and dampers to the
problem. Then, the `Build` method must be called. The methods `SetInitialDisp`
and `SetInitialVel` can then be used to set initial displacements and velocities
to masses. Lastly, the `FixMass` is used to set masses as fixed, i.e. make so
that they have always zero displacement. Finally, the `Integrate` method can be called.

Some post processing methods available are:

- `PrintMassTimeHistory`, which prints to `stdout` the time series of one specific mass's displacement, speed and acceleration. This can, then, be plotted in any csv plotting tool.
- `GetMassMaxAbsAccel` returns max. absolute value of acceleration of a specific mass.
- `GetMassMaxAccel` returns max. value of acceleration of a specific mass.
- `GetMassMinAccel` returns min. value of acceleration of a specific mass.

For more details, see the `software/problem.h` file.

#### Example

```{.cpp caption="Example of using Problem class to perform dynamic simulation. Source: software/problem_test.cc DampedOscillatorPlotTest"}
    Problem p = Problem();
    p.AddMass(1.0, 0.0, 0.0);

    p.AddMass(20.0, 1.0, 0.0);
    p.AddSpring(0, 1, 30.0);
    p.AddDamper(0, 1, 2.9);
    p.Build();
    p.FixMass(0);
    p.SetInitialDisp(1, 1.0);

    p.Integrate(40);

    std::cout << "DampedOscillatorPlotTest output:\n";
    p.PrintMassTimeHistory(1);
```
<div class="largeFig">
![Plot of output of damped oscillator simulation example. Created with chart-studio.plotly.com/](figs/dampedOsc.png){#fig:dampedOsc width=100%}
</div>


## Genetic Algorithm

## Meta-models for mechanical optimization
[@Driemeier_undated-za]
[@Gu2018-uk]
[@Gu2018-tf]
[@Lee2022-uz]
[@Wilt2020-np]

## Physics Informed Neural Networks (PINNs)
[@Nascimento2020-xp]


