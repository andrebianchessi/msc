# Methods {#sec:methods}

## Overview

The method used to achieve the objectives listed at @sec:objectives was to
implement a library that is capable of:

1. Defining arbitrary COPs.

2. Solving COPs with either P-GA or E-GA:
- E-GA evaluates candidate solutions using explicit time integration (see @sec:eti).
- P-GA evaluates candidate solutions using PIMs which describe the position of
each mass as a function of *time* and of the *values of the masses and springs
of the system*.

Using this library, we performed a few COP case studies using both
P-GA and E-GA. The performance and the result of each algorithm were
then compared. *Performance* of each algorithm was measured by the total
processing time needed for the optimization to finish. We compared the *results* by
comparing the maximum acceleration the mass would experience with the optimal
solution found by the algorithm.

## Polynomial Models {#sec:polynomials}

## Physics-Informed Machine Learning Models (PIMs) {#sec:methods_pim}

### Formulation {.unnumbered}

## P-GA {#sec:methods_pga}

## E-GA {#sec:methods_ega}


## Explicit Time Integration Software {#sec:software_eti}

### Implementation {.unnumbered}

The [Problem (~/software/problem.h)](https://github.com/andrebianchessi/msc/blob/7cf80c4f85161acef1c2946262259ad1c3e8f4af/software/problem.h#L17) class, together with other classes it references, encapsulates all the logic related to the dynamic simulation of [CMs](#sec:cms) using ETI.

### Usage {.unnumbered}

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

### Example {.unnumbered}
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

### Implementation {.unnumbered}

[Evolution (~/software/evolution.h)](https://github.com/andrebianchessi/msc/blob/main/software/evolution.h) is a template class that encapsulates the logic related to GA. Note that this template must be of a class that is a child of the [Creature (~/software/creature.h)](https://github.com/andrebianchessi/msc/blob/main/software/creature.h) abstract class (interface).

### Usage {.unnumbered}

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
