# Methods {#sec:methods}

## Overview

The method used to achieve the objectives listed at @sec:objectives was to
implement a library that is capable of:

1. Defining arbitrary COPs.

2. Solving COPs with either [P-GA](#sec:methods_pga) or [E-GA](#sec:methods_ega):
- E-GA evaluates candidate solutions using explicit time integration (see @sec:eti).
- P-GA evaluates candidate solutions using [PIMs](#sec:methods_pim) which describe the position of
each mass as a function of *time* and of the *values of the masses and springs
of the system*. The PIMs used are linear regression models (see @sec:polynomials).

Using this library, we performed a few COP case studies using both
P-GA and E-GA. The performance and the result of each algorithm were
then compared. *Performance* of each algorithm was measured by the total
processing time needed for the optimization to finish. We compared the *results* by
comparing the maximum acceleration the mass would experience with the optimal
solution found by the algorithm.

## Polynomial Models {#sec:polynomials}

The PIMs we used in this work are based on linear regression models.
The expression of the models is defined by the number of springs/dampers of the
[CM](#sec:cms) we're trying to optimized, and by a parameter $h$ -
referred to as the **order** of the models - which defines
the highest order of the monomials.

The models have the following expression:

$$
\begin{aligned}
p_h(t, k_1, k_2, ... , k_i, c_1, c_2,..., c_j) = \sum_{\lambda=0}^{\lambda=h} a_\alpha t^\lambda + &\sum_{Z} a_\alpha t^\beta
k_1^{\gamma_1} k_2^{\gamma_2} ... k_i^{\gamma_i}
c_1^{\omega_1} c_2^{\omega_2} ... c_j^{\omega_j}\\
&Z = \{1 <= \beta < h, (\gamma_i = 0 \text{ OR } \gamma_i = 1), (\omega_j = 0 \text{ OR } \omega_j = 1), {\textstyle(\sum \gamma_i + \sum \omega_j = 1)}\}
\end{aligned}
$$
{#eq:polynomials}

In plain english, that means that the polynomial is a linear combination of:

- All the powers of $t$ from $0$ to $h$
- Cross product of $[t^{h-1}, ..., t^2, t]$ and $[k_1, ..., k_i, c_1, ..., c_j]$

For example, let's consider a [CM](#sec:cms) with 2 springs and 1 damper.
Following are the expression that the models would take:

$$
p_0(t, k_1, k_2, c_1) = a_0
$$

$$
p_1(t, k_1, k_2, c_1) = a_0 + a_1t
$$

$$
p_2(t, k_1, k_2, c_1) = a_0 + a_1t + a_2t^2 + a_3tk_1 + a_4tk_2 + a_5tc_1 
$$

$$
p_3(t, k_1, k_2, c_1) = a_0 + a_1t + a_2t^2 + a_3t^3
+a_4t^2k_1 + a_5t^2k_2 + a_6t^2c_1
+a_7tk_1 + a_8tk_2 + a_9tc_1 
$$

### Reasoning behind the models' architecture

The goal of this project was not to determine optimal model to represent the
dynamic response of [CMs](#sec:cms); but rather just to analyze the
[P-GA](#sec:methods_pga) approach as a whole.

The models defined by @eq:polynomials are convenient for many reasons:

- Their gradients are easy to compute since they're are linear.

- For a given problem, the whole model architecture is defined by a single parameter
(the order) of the models. This makes it super easy to increase/decrease the
model complexity when needed.

- Implementing differentiation and
linear combinations of polynomial models is relatively easy to implement when
compared to other models such as Neural Networks.

Note that the architecture of the models have a significantly strong assumption
behind them: the largest order of **t** is always larger than the order of the other input
variables ($k_i$ and $c_j$). This was characteristic is very much by design,
because we know that the dynamic response of the [CMs](#sec:cms)
are usually not linear with respect to time. Still, the training of the
models should be able to identify which monomials are more significant to model
the systems.

Of course, an immediate opportunity of further research is to use different
models and analyze how that changes the performance of [P-GAs](#sec:methods_pga).

### Automatic differentiation and linear combination {#sec:polynomials_dif}

As is more thoroughly described in @sec:methods_pim, differentiating the models
with respect to time and linearly combining them are two tasks that are required
to build the loss function used to train the models.

As can be seen in [~/software/polynomial.cc](https://github.com/andrebianchessi/msc/blob/main/software/polynomial.cc), the classes implemented for representation/manipulation of
polynomial models support automatic differentiation and the operators of
polynomial instances [Poly and Polys (~/software/polynomial.cc)](https://github.com/andrebianchessi/msc/blob/main/software/polynomial.cc) have been implemented so that they can be linearly
combined even through matrix multiplications (see [this test](https://github.com/andrebianchessi/msc/blob/897a324b2d0e1f0b12d1c211f10f0cee64fd2f7c/software/polynomial_test.cc#L400) for example). 

#### Side note {.unnumbered}

On a side note, I'd recommend the reader to skim through the implementation of
[~/software/polynomial.cc](https://github.com/andrebianchessi/msc/blob/main/software/polynomial.cc),
as implementing these classes was one of the intermediary tasks of this project
which I found the most challenging and interesting.

### Software

[~/software/polynomial.h](https://github.com/andrebianchessi/msc/blob/main/software/polynomial.h)
contains all the code related to the polynomial models. The [Poly (~/software/polynomial.h)](https://github.com/andrebianchessi/msc/blob/main/software/polynomial.h) class defines a single instance
of a polynomial, and the [Polys (~/software/polynomial.h)](https://github.com/andrebianchessi/msc/blob/main/software/polynomial.h) class handles the linear combination of *Poly* instances.

As usual, [~/software/polynomial_test.cc](https://github.com/andrebianchessi/msc/blob/main/software/polynomial_test.cc) contains tests which serve as documentation and usage examples.


## Physics-Informed Machine Learning Models (PIMs) {#sec:methods_pim}

### Architecture of the models {#sec:methods_pim_models}

As seen in @sec:pim, *Physics-Informed Machine Learning Models* are trained so that they approximate the solution
to a differential equation. In this work, we're dealing not with **one**, but with a **system** of differential equations that describe the [COP](#sec:cop).

Consider an arbitrary [COP](#sec:cop) of $n$ masses. As seen in @sec:eti, the $n$ equations that describe the system, as obtained with the [Discrete element method](@sec:dem),
have the following form:

$$
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
{#eq:methods-pim-ode}

To use *Physics-Informed Machine Learning Models* for this problem, our approach was to have **one model per mass** that approximates
the **displacement ($x_n$) of each mass as a function of time ($t$) and of the constants of the springs and dampers ($k_0, \cdots,k_i, c_0, \cdots,c_j$)**. I.e. for a system of $n$
masses, $i$ springs and $j$ dampers, we defined the [polynomial models](#sec:polynomials) $P_0(t, k_0, \cdots,k_i, c_0, \cdots,c_j), \cdots, P_n(t, k_0, \cdots,k_i, c_0, \cdots,c_j)$ so that $P_0$ models $x_0$, $P_1$ models $x_1$, etc.

Note that the models can easily be differentiated
with respect to time so that we can obtain the velocity and acceleration of each mass.
The order of the polynomials ($h$) was a hyperparameter chosen for the experiments
(see @sec:experiments for more details).

### Time discretization {#sec:methods_pim_t_disc}

Given that the time responses of the [CMs](@sec:cms) are usually very non-linear with
respect to time (see @fig:dampedOsc for example), using a model that describes
the displacement of a mass for the whole
duration of the impact would not be efficient. The model would need to have a very
high order, which can make the training very slow.

To solve that, the approach we took was to discretize the time into multiple "buckets".
Let's say we want models that describe a [COP](#sec:cop) from $t=0$ to $t=T$. Instead of
having a set of models (one for each mass) that describes the the displacement of each mass
as a function of time from $t=0$ to $t=T$, we created a set of models (one for each mass) that describe
the displacement of the masses from $t=0$ to $t=T_0$, then another set of models for $t=T_0$ to $t=T_1$,
and so on until a set of models for $t=T_i$ to $t=T$. The number of "time buckets" was another
hyperparameter chosen for the experiments
(see @sec:experiments for more details).

To train those models, we start by training the first set of models - they describe the displacement
of the masses from $t=0$ to $t=T_0$. Let's call these the $t_0$ models. The initial conditions
(the displacement and velocity of each mass) are given by the
[COP](#sec:cop) statement. Then, to train the next set of models - the $t_1$ models - we used the $t_0$ models to find
the conditions (displacement and velocity of each mass) at $t=T_0$. These conditions are considered
initial conditions for the next set of models. This process continues until all the $t_i$ models are trained. At the end of this process,
we have a set of models for each "time bucket". See [Pimodels::Train (~/software/pimodel.cc)](https://github.com/andrebianchessi/msc/blob/e24929a15dad217a5f2366a9544765ea25eea6f5/software/pimodel.cc#L616) for more detail.

Note that the models have as input not only the time, but also the values of the springs and dampers.
When using a "previous" set of models to determine the initial conditions to train the "next" set of models,
we choose the intermediary values for the springs and dampers. I.e. for every spring/damper that can have
values from $a$ to $b$, we used $(a+b)/2$ to evaluate the models.

### Normalization {#sec:methods_pim_normalization}

Normalizing the inputs to the models can yield faster training [@MlBook, p. 365].
Thus, the input to all models were the normalized values from `0` to `1`.
For the springs, `0` corresponded to minimum possible value of the springs elastic constant and `1` to the maximum.
The analogous was used for the dampers. For time, `0` corresponded to the start of the impact and `1` to the
end of the impact (`t=T`). A linear normalization was used.


### Training: Minimizing the loss function {#sec:methods_pim_training}

As described in @sec:methods_pim_t_disc, sets of models are progressively trained for each "time bucket".
This section explains how each set of models is trained.

Initially, all the models are created with all the coefficients equal to zero; i.e. all the
polynomial coefficients are 0.

*Training* is the name of the process used to find optimal values for the parameters of the models. It is done by minimizing, through Stochastic Gradient Descent [@MlBook, p. 184], a loss function.

Following the usual formulation of Physics Informed Machine learning [@Thuerey2021-ut],
the loss function is composed of 3 parts:

- $L_{x}$: Initial displacement loss

- $L_{\dot{x}}$: Initial velocity loss

- $L_{\ddot{x}}$: Physics loss

$L_{x}$ measures how well the models estimate the initial displacement of the masses.
$L_{\dot{x}}$ measures how well the models estimate the initial velocity of the masses.
$L_{\ddot{x}}$ measures how well the models obey the physics, i.e. how well the fit the [COP](#sec:cop)'s ODE (@eq:methods-pim-ode).

For the following text, consider that $s$:

$$
s={k_0}_{norm}, ...,{k_i}_{norm}, {c_0}_{norm}, ...,{c_j}_{norm}
$$
{#eq:s}

represents a set of springs and dampers that define a solution to a [COP](#sec:cop). The $_{norm}$ subscript indicates that the values are normalized considering the maximum and
minimum possible values of each spring and damper (see @sec:methods_pim_normalization).
As an example, let's consider a
[COP](#sec:cop) comprized of two masses connected by a spring of elastic constant $k_0$ and a damper with damping coefficient $c_0$. The spring can have values from $10$ to $20$,
and the damper can have values from $50$ to $100$. $s=\{0, 0.5\}$ represents a system with $k_0 = 10$ and $c_0 = 75$.

In our chosen approach, we wan the models to not only be a function of time, but also of the system's
parameters (the springs and the dampers). Thus, we want the minimization of the losses to force the models
to work well for multiple values of springs and dampers.
We achieved that by computing the losses for multiple possible values of the springs and dampers as follows:

Given the hyperparameter $\alpha$, we first define a random set $S$ of size $\alpha$ that contains random values of $s$ (see @eq:s):

$$
S = \{s_1, s_2, ... , s_{\alpha} \}
$$

For a [COP](#sec:cop) of $n$ masses, let $P_i$ indicate the model of the mass $i$ (see @sec:methods_pim_models), $L_{x}$ is then defined as:

$$
L_{x} = \sum_{i=1}^{n}\sum_{j=1}^{\alpha} (P_i(t = 0, s_j) - x_i|_{t=0})^2
$$
{#eq:lx}

In plain english, that equation is read as: Sum of the squared error of the model's predicted displacement and the actual displacement at $t=0$, for all the masses
using all the random values of springs and dampers we randomly picked.

Note that the $x_i|_{t=0}$ are the initial conditions of each mass, which are part of the [COP](#sec:cop) definition.

The definition of $L_{\dot{x}}$ is analogous, but for the initial velocities:

$$
L_{\dot{x}} = \sum_{i=1}^{n}\sum_{j=1}^{\alpha} (\dot{P_i}(t = 0, s_j) - \dot{x_i}|_{t=0})^2
$$
{#eq:ldotx}

@eq:methods-pim-ode is the ODE that the $\ddot{x}_i$ must satisfy. This equation provides all the accelerations as a function of the springs, dampers, and of the displacement and velocities of all masses:

$$
\ddot{x_i} = F_i(x_0, x_1, \cdots, x_n, \dot{x_0}, \dot{x_1}, \cdots, \dot{x_i}, k_0, k_1, \cdots,c_0, c_1, \cdots)
$$
{#eq:xidotdot}

We want the models to fit @eq:xidotdot not only for $t=0$, but for all the interval of the [COP](#sec:cop).
Consider that $s_t$:

$$
s_t = t_{norm}, {k_0}_{norm}, ...,{k_i}_{norm}, {c_0}_{norm}, ...,{c_j}_{norm}
$$
{#eq:st}

represents a set of springs and dampers that define a solution to a [COP](#sec:cop) and a time from 0 to 1.
This is basically the same definition of @eq:s, but with *time* also as an argument. Similarly to how we built the other
losses, given the hyperparameter $\beta$,
we first define a random set $S_{t}$ of size $\beta$ that contains random values of $s_t$ (see @eq:st):

$$
S_{t} = \{{s_t}_1, ... , {s_t}_{\beta} \}
$$

The physics loss function is then defined as:

$$
L_{\ddot{x}} = \sum_{i=1}^{n}\sum_{j=1}^{\beta} (\ddot{P_i}({s_t}_j) -
F_i(
    P_0({s_t}_j), P_1({s_t}_j), \cdots,
    \dot{P_0}({s_t}_j), \dot{P_1}({s_t}_j), \cdots,
    {s_t}_j
    ))^2
$$
{#eq:ldotdotx}

That equation is basically the sum of the errors of the second derivative of the model
and what that second derivative *should* be according to the system's ODE.
Since the ODE expects the displacements and velocities of all the masses, we use the models themselves
as estimators for those.

The total loss that we want to minimize is simply the sum of all losses:

$$
L = L_x + L_{\dot{x}} + L_{\ddot{x}}
$$
{#eq:loss}

#### Change of variables for the time derivatives {.unnumbered}

The first element of the residue from @eq:ldotdotx - $\ddot{P_i}({s_t}_j)$ - is a
second derivative of the model with respect to time. As explained in @sec:methods_pim_normalization,
the input of the model is not *the actual time* but rather the normalized time (from 0 to 1).
To clearly distinguish those two, in this section $T_g$ represents the *global* time, and $T_l$ represents
the *local* time, which is normalized from $0$ to $1$. Using that notation, @eq:ldotdotx becomes:

$$
L_{\ddot{x}} = \sum_{i=1}^{n}\sum_{j=1}^{\beta} (\frac{d^2P_i}{dT_g^2}({s_t}_j) -
F_i(
    P_0({s_t}_j), P_1({s_t}_j), \cdots,
    \frac{dP_0}{dT_g}({s_t}_j), \frac{dP_1}{dT_g}({s_t}_j), \cdots,
    {s_t}_j
    ))^2
$$

For more clarity, let's turn our attention to the residue that is being summed:
$$
\frac{d^2P_i}{dT_g^2}({s_t}_j) -
F_i(
    P_0({s_t}_j), P_1({s_t}_j), \cdots,
    \frac{dP_0}{dT_g}({s_t}_j), \frac{dP_1}{dT_g}({s_t}_j), \cdots,
    {s_t}_j
    )
$$

With the chain rule, that becomes:
$$
\left( \frac{dT_l}{dT_g} \right)^2\frac{d^2P_i}{dT_l^2}({s_t}_j) -
F_i(
    P_0({s_t}_j), P_1({s_t}_j), \cdots,
    \frac{dT_l}{dT_g}\frac{dP_0}{dT_l}({s_t}_j), \frac{dT_l}{dT_g}\frac{dP_1}{dT_l}({s_t}_j), \cdots,
    {s_t}_j
    )
$$

As we use smaller "time buckets", the $dT_l/dT_g$ term increases. For example, let's say we're
considering a total impact duration of 0.05 seconds, and using a *time discretization* of 10
(i.e. 10 time buckets). In this case we'll first train a set of models for $T_g=0$ to $T_g=0.005$;
then another set of models from $T_g=0.005$ to $T_g=0.010$ and so on.
For the first set of models, $T_l = 200*T_g$, so $dT_l/dT_g = 200$. That derivative is the same for all
other set of models.

Thus, it's easy to see that $(dT_l/dT_g)^2$ rapidly increases as we discretize the time.
This causes the values and gradients of $L_{\ddot{x}}$ to "blow up", which makes the training process take longer.

Therefore, it's convenient to rewrite the residue as follows:

$$
\left( \frac{dT_l}{dT_g} \right)^2 \left( \frac{d^2P_i}{dT_l^2}({s_t}_j)
- \left(\frac{dT_g}{dT_l} \right)^2
F_i(
    P_0({s_t}_j), P_1({s_t}_j), \cdots,
    \frac{dT_l}{dT_g}\frac{dP_0}{dT_l}({s_t}_j), \frac{dT_l}{dT_g}\frac{dP_1}{dT_l}({s_t}_j), \cdots,
    {s_t}_j
    )
\right)
$$

In Stochastic Gradient Descent, only one of the loss terms is evaluated at a time; and the gradient of that
loss term is used to update the parameter of the models. Since that gradient is scaled anyway, we can
disconsider the $\left(dT_l/dT_g \right)^2$ term when computing the gradient.

This alternative formulation contains the $\left(dT_g/dT_l \right)^2$, which will decrease rapidly
as we use higher *time discretization*. Still, in the experiments we analyzed we saw that this formulation
behaves much better numerically than the former one.

To summarize, the physics loss can be expressed as:

$$
L_{\ddot{x}} = \sum_{i=1}^{n}\sum_{j=1}^{\beta}
\frac{d^2P_i}{dT_l^2}({s_t}_j)
- \left(\frac{dT_g}{dT_l} \right)^2
F_i(
    P_0({s_t}_j), P_1({s_t}_j), \cdots,
    \frac{dT_l}{dT_g}\frac{dP_0}{dT_l}({s_t}_j), \frac{dT_l}{dT_g}\frac{dP_1}{dT_l}({s_t}_j), \cdots,
    {s_t}_j
    )
$$
{#eq:ldotdotxDiff}

The $dP_i/dT_l$ derivatives are easy to compute because they're simple differentiations of polynomials,
and the $dT_l/dT_g$ is also trivially computed simply with the expression that normalizes the time.
The $F_i$ is basically a linear combination of models. For these reasons, the automatic differentiation
and linear combination of models was necessary (see @sec:polynomials_dif).

### Example: Putting it all together {#sec:methods_pim_example}

### Software {#sec:methods_pim_software}

## P-GA {#sec:methods_pga}

## E-GA {#sec:methods_ega}


## Explicit Time Integration Software {#sec:software_eti}

### Implementation

The [Problem (~/software/problem.h)](https://github.com/andrebianchessi/msc/blob/7cf80c4f85161acef1c2946262259ad1c3e8f4af/software/problem.h#L17) class, together with other classes it references, encapsulates all the logic related to the dynamic simulation of [CMs](#sec:cms) using ETI.

### Usage

[~/software/problem_test.cc](https://github.com/andrebianchessi/msc/blob/main/software/problem_test.cc) contains many examples of how the software we
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

For more details, see [~/software/problem.h](https://github.com/andrebianchessi/msc/blob/main/software/problem.h).

### Example
[DampedOscillatorPlotTest (~/software/problem_test.cc)](https://github.com/andrebianchessi/msc/blob/e7e048d554f82161702b1f90b3878957dbb0538b/software/problem_test.cc#L684)

```{.cpp caption="Example of using Problem class to perform dynamic simulation."}
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

![Plot of output of damped oscillator simulation example. Created with gnuplot](figs/dampedOsc.png){#fig:dampedOsc width=100% style="scale:1;"}

## Genetic Algorithm Software {#sec:software_ga}

### Implementation

[Evolution (~/software/evolution.h)](https://github.com/andrebianchessi/msc/blob/main/software/evolution.h) is a template class that encapsulates the logic related to GA. Note that this template must be of a class that is a child of the [Creature (~/software/creature.h)](https://github.com/andrebianchessi/msc/blob/main/software/creature.h) abstract class (interface).

### Usage

#### Arbitrary problem {.unnumbered}
To perform an optimization using GA, the first step is to define a child class of
[Creature (~/software/creature.h)](https://github.com/andrebianchessi/msc/blob/main/software/creature.h).
The only definitions required are the ```dna``` attribute and the ```GetCost``` function.
[C (~/software/creature_test.h)](https://github.com/andrebianchessi/msc/blob/3820ac9bd39e60e2aed403fba2227cd116772228/software/creature_test.h#L18) class is an example. The creature class, in this example, represents a candidate $\{x,y\}$ pair that minimizes $f(x,y) = x^2 + y^2 + 2x + y$.

Once the *child creature* class is defined, an `Evolution` object can be instantiated and used to search for optimal solutions. [EvolveTest (~/software/evolution_test.cc)](https://github.com/andrebianchessi/msc/blob/3820ac9bd39e60e2aed403fba2227cd116772228/software/evolution_test.cc#L326) contains an example of how that's done.

#### COP {.unnumbered}

The *child class* for [CMs](#sec:cms) is already defined at [ProblemCreature (~/software/problem_creature.h)](https://github.com/andrebianchessi/msc/blob/main/software/problem_creature.h). Its constructor uses an auxiliary class [ProblemDescription (~/software/problem_description.h)](https://github.com/andrebianchessi/msc/blob/main/software/problem_description.h), which allows us to easily describe a COP as described in @sec:cop. The extra two parameters from its constructor specify the mass whose maximum acceleration we want to minimize and the length of the simulation we'll perform. The creature's loss function is the maximum absolute acceleration that mass will suffer during the dynamic simulation.

#### Example {.unnumbered}

[EvolutionUntilConvergenceTest (~/software/problem_creature_test.cc)](https://github.com/andrebianchessi/msc/blob/d69d36a973dce9807674b023ad8bfd05b2b7a612/software/problem_creature_test.cc#L89) contains an example in which we find values for the springs and dampers of system at @fig:crashworthiness that minimize the maximum acceleration that $m_5$ would suffer if the system was moving with a constant speed from right to left and hit an immovable wall on the left. A simplified version of the code is listed below. Note that the parameters we pass to the ```Evolve``` method determine our stop condition and if the results should be printed to ```stdout```. Some values of the best solution found are listed after it, and @fig:msdsGa shows how the sum of the loss of the fittest population progresses with the generations.

```{.cpp caption="Example of how to use the code we wrote to solve COPs from @fig:crashworthiness"}
ProblemDescription pd = ProblemDescription();
pd.AddMass(1.0, 0.0, 0.0);  // m0
pd.AddMass(300, 1.0, 1.0);  // m1
pd.AddMass(120, 1.0, 0.0);  // m2
pd.AddMass(150, 1.0, 3.0);  // m3
pd.AddMass(700, 2.0, 0.0);  // m4
pd.AddMass(80, 3.0, 0.0);   // m5

double min = 100.0;
double max = 100000;
pd.AddSpring(0, 1, min, max);
pd.AddSpring(1, 2, min, max);
pd.AddSpring(1, 3, min, max);
pd.AddSpring(1, 4, min, max);
pd.AddDamper(1, 4, min, max);
pd.AddSpring(0, 2, min, max);
pd.AddDamper(0, 2, min, max);
pd.AddSpring(2, 4, min, max);
pd.AddDamper(2, 4, min, max);
pd.AddSpring(0, 3, min, max);
pd.AddDamper(0, 3, min, max);
pd.AddSpring(3, 4, min, max);
pd.AddDamper(3, 4, min, max);
pd.AddSpring(4, 5, min, max);
pd.AddDamper(4, 5, min, max);

pd.SetFixedMass(0);
pd.AddInitialVel(200.0);

// Create population of 20 creatures
std::vector<ProblemCreature> pop = std::vector<ProblemCreature>();
for (int i = 0; i < 20; i++) {
    pop.push_back(ProblemCreature(&pd, 5, 0.15));
}

// Find optimal solutions
Evolution<ProblemCreature> evolution = Evolution<ProblemCreature>(&pop);
double cost0 = evolution.FittestCost();
auto p = evolution.Evolve(0.01, true);

// Check values of best solution
Problem best = pd.BuildFromDNA(evolution.GetCreature(0)->dna).val;
print("k1: ", best.springs[0].Get_k());
print("k2: ", best.springs[1].Get_k());
print("k3: ", best.springs[2].Get_k());
print("k4: ", best.springs[3].Get_k());
print("c4: ", best.dampers[0].Get_c());
print("k5: ", best.springs[4].Get_k());
print("c5: ", best.dampers[1].Get_c());
```

Some of the values of springs and dampers of the optimal solution we found are:

| Component|   Value  |
|:--------:|:--------:|
|        k1|   28870.3|
|        k2|     73618|
|        k3|   43417.4|
|        k4|   70919.7|
|        c4|    4957.4|
|        k5|   21287.2|
|        c5|    4957.4|

![Example solution of COP: Sum of loss function of fittest population vs Generation](figs/msdsGa.png){#fig:msdsGa width=80% style="scale:1;"}

## Experiments {#sec:experiments}
