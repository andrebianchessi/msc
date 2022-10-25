## Genetic Algorithm

*The genetic algorithm is an optimization and search technique based on the principles of genetics and natural selection* [@Haupt2004-hj]. Some of its advantages, when compared to gradient base methods, are [@Haupt2004-hj]:

- Doesn’t require derivative information
- Simultaneously searches from a wide sampling of the cost surface
- Deals with a large number of variables
- Is well suited for parallel computers
- Optimizes variables with extremely complex cost surfaces (they can jump out of local minima)

There's no single recipe for success of GAs. Many details vary slightly between
implementations - such as the use of *elitism* the criteria for selecting mates
and the breeding process - but @fig:ga illustrates a high level overview of all
components of the Genetic Algorithm as we implemented it. In this section, the
components, as we implemented them in our code, are briefly explained and an
example involving all the steps is presented.

The main sources this section is based on are [@Haupt2004-hj; @Lam2021-gp;
@Lam2021-hg].

![GA flowchart. Source: Author](figs/ga.png){#fig:ga scale=1 style="scale:1;"}

The [Evolve (~/software/evolution.tcc)](https://github.com/andrebianchessi/msc/blob/e7e048d554f82161702b1f90b3878957dbb0538b/software/evolution.tcc#L203) method is the top-level-method which executes the GA optimization. It basically performs the initial sorting of the population, checks for convergence and successively calls the [step (~/software/evolution.tcc)](https://github.com/andrebianchessi/msc/blob/e7e048d554f82161702b1f90b3878957dbb0538b/software/evolution.tcc#L175) method, which performs one iteration of optimization of the population.

### GA steps

#### Choose DNA that represents creature {.unnumbered}

Given GA's biological inspirations, possible solutions to the problem at hand are often referred to as *creatures*. For example, if we use GA to find the minimum value of a function $f(x,y)$, 
$\{x,y\}$ pairs, such as $\{0,0\}$, $\{1,0\}$, can be thought of as *creatures*.

In GA, we must define what the DNA of our problem's creatures is.
For *real-encoded* - also known as *continuous* - GA, which is the one we're interested in, the DNA is simply a vector of real numbers. So for the problem we mentioned above in which creatures are described by pairs of $x$ and $y$, we could choose that the creatures' DNA is simply a vector in which the first element is the value of $x$ and the second is the value of $y$. Note that we could also choose the other order. Thus, it's up to the user to choose how to model the problem.

Another example: let's say the problem at hand is that of finding values for the springs and dampers of the system illustrated at @fig:discreteElementSimple4 which minimize the maximum acceleration $m_3$ suffers if the whole system is traveling with a constant speed towards the left and hits an immovable object.
A *creature* in this problem is a set of values of springs' and dampers' constants.
We can choose, for example, that the DNA that represents a creature is the vector $[k_{01},c_{01},k_{12},c_{12},k_{13},c_{13},k_{03},c_{03}]$.

Note that when defining a DNA, we must also define what our domain is, i.e. in what range each DNA member must be in.

#### Choose loss function {.unnumbered}

In this step, we need to define a function we want to minimize. This function takes a creature - or more specifically its DNA - as input, and returns a real value, i.e. a scalar.

#### Choose hyperparameters {.unnumbered}

Hyperparameters are top-level parameters we choose for the optimization algorithm itself. In this case, we must choose:

- Population size
- Population survival rate
- Mutation rate

The meaning of these values will become clearer in the following example.

#### Create initial population {.unnumbered}

At this step, we initialize a population of random creatures.

#### Sort population and remove the less fit {.unnumbered}

We begin this step by calculating the loss function of each creature. Then, we sort the population by ascending loss function value. A portion, determined by the *population survival rate* hyperparameter, of the less fit population (the creatures with the highest loss function values) is then removed from the population. This stage is equivalent to the survival of the fittest in the evolutionary process.

#### Select mates {.unnumbered}

In this step, creatures which will mate to create offspring that will replace the creatures which were removed at the last stage must be selected.

There are many different approaches to doing this. Some of them are highlighted at [@Lam2021-gp;
@Lam2021-hg]. In our code, we used what's known as *Biased Roulette Wheel Selection*. In this process, two creatures are selected randomly, but the probability of a creature being selected is proportional to how low their loss function value is. For more detail on how this done, we encourage readers to look at the implementation of [getParents (~/software/evolution.tcc)](https://github.com/andrebianchessi/msc/blob/e7e048d554f82161702b1f90b3878957dbb0538b/software/evolution.tcc#L110).

#### Mate {.unnumbered}

In this step, two creatures which were selected in the previous step - the parents - mate to create two new creatures : the children. Again, there are also many different techniques for this process (see @Haupt2004-hj). We implemented what's known as *Radcliff blending method* @Haupt2004-hj, which is as follows:

Let $\text{D}_{p0}$ and $\text{D}_{p1}$ be the DNA of the parents, and $\text{D}_{c0}$ and $\text{D}_{c1}$ be the DNA of the children which will be created.
$\text{D}_{x}[i]$ corresponds to the $i$-th position of $x$'s DNA.

For every valid index $i$:

1. A random variable $\beta$ in the range $[0, 1]$ is created.
2. $\text{D}_{c0}[i] = \beta \text{D}_{p0}[i] + (1-\beta \text{D}_{p1}[i])$
3. $\text{D}_{c1}[i] = \beta \text{D}_{p1}[i] + (1-\beta \text{D}_{p0}[i])$

This logic is implemented at the [Mate (~/software/creature.cc)](https://github.com/andrebianchessi/msc/blob/e7e048d554f82161702b1f90b3878957dbb0538b/software/creature.cc#L19) method.

The two created children fill the gaps in the population of the less fit creatures which were removed at the *Sort population and remove the less fit* step. If there's only room left in the population for one creature, one of them is just discarded.

#### Done mating? {.unnumbered}

If the current population size is still smaller than the *population size* hyperparameter, we continue with the mating process.

#### Mutate {.unnumbered}

At this stage, we apply random mutations to only the children which were created (which is called *elitism*). *Elitism* his guarantees that the fittest creatures from one iteration of the algorithm are either just as fit or less fit than the fittest ones from the next iteration of the algorithm.

We implemented uniform random mutation, which means that to perform a mutation we do the following:

1. Pick a random child which was created
2. Choose a random position of the DNA
3. Let $a$ and $b$ be the minimum and maximum value acceptable for this position at the DNA. Replace the value at that position with a random number in the interval $[a,b]$

The number of mutations we perform is determined by the *mutation rate* hyperparameter. For a mutation rate $m$, population size $p$, *survival rate* $s$ and the creatures' DNA size $d$, the number of mutations that will be performed is given by $m \cdot p(1-s)d$.

The idea behind this is that we allow us only to mutate the children. Thus, there are $p(1-s)$ creatures that can suffer mutation. Each has $d$ DNA slots. Hence, there are $p(1-s)d$ DNA slots we can mutate. We multiply that by the *mutation rate*, which is a number between $0$ and $1$, and get the number of DNA slots we'll mutate.

This logic is implemented at the [mutate (~/software/evolution.tcc)](https://github.com/andrebianchessi/msc/blob/e7e048d554f82161702b1f90b3878957dbb0538b/software/evolution.tcc#L140) method.

#### Continue? {.unnumbered}

Whenever we reach this stage, we consider to have finished *a generation*. We
now choose to stop if the algorithm has converged, i.e. the value of the loss
function of the fittest creatures is practically the same it was on the previous
generation, or if too many generations have been tried but the algorithm still
hasn't converged.

### GA steps illustrative example

Suppose we want to find the value of $\{x,y\}$ in the region $x \in [-2,2]$ and $y \in [-2,2]$ that minimizes the function

$$
f(x,y) = x^2 + y^2 + 2x + y
$$

#### Choose DNA that represents creature {.unnumbered}

The creatures of this problem are pairs $\{x,y\}$, so we can choose the DNA of each creature to be the vector $[x,y]$

#### Choose loss function {.unnumbered}

Since we want to minimize $f$, the loss function for a creature with DNA $[x_c,y_c]$ can simply be:

$$
f(x_c,y_c) = x_c^2 + y_c^2 + 2x_c + y_c
$$

#### Choose hyperparameters {.unnumbered}

Since this is just an illustrative example, we chose:

- Population size $=4$
- Population survival rate $=0.5=50\%$
- Mutation rate $=0.25$

Note that outside illustrative examples the population sizes are usually much
higher and the mutation rates are usually much smaller.

#### Create initial population {.unnumbered}

By picking random numbers between $-2$ and $2$ we obtained the following population:

$$
[[0.4, 0.16], [-0.9, 0.2], [-0.2, 0.2], [0.1, -0.1]]
$$
{#eq:gaInitialRandomPop}

#### Sort population and remove the less fit {.unnumbered}

The loss function calculated for each creature at @eq:gaInitialRandomPop has the following values:

$$
[1.1456,-0.75,-0.12,0.12]
$$

Thus, the sorted population becomes:
$$
[[-0.9, 0.2], [-0.2, 0.2], [0.1, -0.1], [0.4, 0.16]]
$$

With a survival rate of $50\%$, we get:
$$
[[-0.9, 0.2], [-0.2, 0.2]]
$$

#### Select mates {.unnumbered}

In this case, since there are only 2 creatures, we have no choice but to select them as parents. However, the roulette wheel algorithm, which is implemented at the [getParents (~/software/evolution.tcc)](https://github.com/andrebianchessi/msc/blob/e7e048d554f82161702b1f90b3878957dbb0538b/software/evolution.tcc#L110) method would be as follows:

First we transform the loss value into a fitness value. We can do that by adding $1.75$ to all the loss values, and then taking the inverse of the value:

$$
[-0.75,-0.12] \rightarrow [1,1.63] \rightarrow [1,0.61]
$$

We then normalize the values:
$$
[1,0.61] \rightarrow [0.62, 0.38]
$$

Lastly, we accumulate the values, so that they're always growing.
$$
[0.62, 0.38] \rightarrow [0.62, 1]
$$

Now, to select a parent we create a random number between $0$ and $1$. If the random number is $<=0.62$, we select the first creature as parent. Else, we select the second. We then pick another random number to choose the other parent.
In cases when the same parent is selected twice by chance, we take the next parent.

#### Mate {.unnumbered}

For the random number $\beta = 0.1$, the DNAs of the children would be:

$$
\begin{aligned}
&[0.1\cdot(-0.9) + (1-0.1)\cdot(-0.2), 0.1\cdot(0.2) + (1-0.1)\cdot(0.2)] = [-0.27, 0.2]\\
&[0.1\cdot(-0.2) + (1-0.1)\cdot(-0.9), 0.1\cdot(0.2) + (1-0.1)\cdot(0.2)] = [-0.83, 0.2]
\end{aligned}
$$

After adding them to the population we get:

$$
[[-0.9, 0.2], [-0.2, 0.2], [-0.27, 0.2], [-0.83, 0.2]]
$$


#### Mutate {.unnumbered}

Two children were created. They each have 2 DNA slots. Thus, there are, in total, 4 DNA slots which can suffer mutation. A *mutation rate* of 25\% means we'll mutate 1 DNA slot.

To perform the mutations we first pick a random child. Then, we pick a random DNA index, and replace the value there with a random one between $-2$ and $2$.

For example, with the following events:

1. First child was randomly selected
2. Second DNA index was randomly selected
3. 0.7 was randomly selected in the range [-2,2]

The final population becomes:

$$
[[-0.9, 0.2], [-0.2, 0.2], [-0.27, 0.7], [-0.83, 0.2]]
$$

### Software

#### Implementation {.unnumbered}

[Evolution (~/software/evolution.h)](https://github.com/andrebianchessi/msc/blob/main/software/evolution.h) is a template class that encapsulates the logic related to GA. Note that this template must be of a class that is a child of the [Creature (~/software/creature.h)](https://github.com/andrebianchessi/msc/blob/main/software/creature.h) abstract class (interface).

#### Usage {.unnumbered}

##### Arbitrary problem {.unnumbered}
To perform an optimization using GA, the first step is to define a child class of
[Creature (~/software/creature.h)](https://github.com/andrebianchessi/msc/blob/main/software/creature.h).
The only definitions required are the ```dna``` attribute and the ```GetCost``` function.
[C (~/software/creature_test.h)](https://github.com/andrebianchessi/msc/blob/3820ac9bd39e60e2aed403fba2227cd116772228/software/creature_test.h#L18) class is an example. The creature class, in this example, represents a candidate $\{x,y\}$ pair that minimizes $f(x,y) = x^2 + y^2 + 2x + y$.

Once the *child creature* class is defined, an `Evolution` object can be instantiated and used to search for optimal solutions. [EvolveTest (~/software/evolution_test.cc)](https://github.com/andrebianchessi/msc/blob/3820ac9bd39e60e2aed403fba2227cd116772228/software/evolution_test.cc#L326) contains an example of how that's done.

##### MSDS optimization {.unnumbered}

The *child class* for MSDSs is already defined at [ProblemCreature (~/software/problem_creature.h)](https://github.com/andrebianchessi/msc/blob/main/software/problem_creature.h). Its constructor uses an auxiliary class [ProblemDescription (~/software/problem_description.h)](https://github.com/andrebianchessi/msc/blob/main/software/problem_description.h), which allows us to easily describe an MSDS optimization problem as described in @sec:problemsWeSolved. The extra two parameter from its constructor specify the mass whose maximum acceleration we want to minimize and the length of the simulation we'll perform. The creature's loss function is the maximum absolute acceleration that mass will suffer during the dynamic simulation.

[EvolutionUntilConvergenceTest (~/software/problem_creature_test.cc)](https://github.com/andrebianchessi/msc/blob/d69d36a973dce9807674b023ad8bfd05b2b7a612/software/problem_creature_test.cc#L89) contains an example in which we find values for the springs and dampers of system at @fig:crashworthiness that minimize the maximum acceleration that $m_5$ would suffer if the system was moving with a constant speed from right to left and hit an immovable wall on the left. A simplified version of the code is listed below. Note that the parameters we pass to the ```Evolve``` method determine our stop condition and if the results should be printed to ```stdout```. Some values of the best solution found are listed after it, and @fig:msdsGa shows how the sum of the loss of the fittest population progresses with the generations.

```{.cpp caption="Example of how to use the code we wrote to optimize MSDS from @fig:crashworthiness"}
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

![Example optimization of MSDS: Sum of loss function of fittest population vs Generation](figs/msdsGa.png){#fig:msdsGa width=80% style="scale:1;"}