# Literature Review

## Crashworthiness models (CMs) {#sec:cms}

Crashworthiness models (CMs) are models of vehicles used to analyze the safety of
their occupants in a crash. @Fig:crashworthiness shows an example. They can be
very elaborate, but very simple ones are also used on initial phases of vehicle
design. These simple ones are one-dimensional Masses-Springs-Dampers Systems,
i.e. they're one-dimensional and are comprised of only ideal masses,
springs and dampers. They can be used to estimate the maximum acceleration a
driver would suffer in a crash. Alternatively, we can look for values of the
springs and dampers that would minimize the maximum acceleration the drive would
suffer in a crash.

![Example of Crashworthiness model. @Tbl:crashworthiness has the legend. Source: [@Mostafa2011-kc]](figs/crashworthiness.png){#fig:crashworthiness width=80% style="scale:1;"}


|  Mass No. |                    Vehicle components |
|:---------:|--------------------------------------:|
|          1|                    Engine and Radiator|
|          2|             Suspension and Front Rails|
|          3|             Engine Cradle and Shotguns|
|          4| Fire Wall and Part of Body on Its Back|
|          5|                               Occupant|
: Legend of @Fig:crashworthiness. Source: [@Mostafa2011-kc] {#tbl:crashworthiness}

## Crashworthiness optimization problem (COP) {#sec:cop}

A *Crashworthiness optimization problem* (COP) is the optimization problem of a
[CM](#sec:cms).

The main objective of this study was to analyze the cost/benefits of
using *P-GAs* vs *E-GAs* (see @sec:objectives) for optimization
of mechanical systems. To achieve this objective, we analyzed multiple
COPs using each algorithm and compared
the results. The following section defines the problem statement of a COP:

### Problem Statement

Consider a mechanical system comprized of ideal masses,
ideal linear springs and ideal linear dampers, such as from @fig:copSystem.
($m_0$, $m_1$, ...,$m_n$) represent the masses,
($k_0$, $k_1$, ...,$k_i$) represent the elastic constants of the springs and
($c_0$, $c_1$, ...,$c_j$) represent the damping coefficient of the dampers.
Note that $m_0$ is fixed, but all the others have **arbitrary initial displacement
and velocities**.

The optimization problem is stated as:
**Given the masses, the initial conditions ($x_1(t=0)$, ...,$x_n(t=0)$, $\dot{x_1}(t=0)$, ...,$\dot{x_n}(t=0)$),
the maximum and minimum values of each $k_i$ and $c_j$,
and an impact duration $T$, find ($k_0$, ...,$k_i$, $c_0$, ...,$c_j$) that minimize the maximum $\ddot{x_n}$ that will occur from $t=0$ to $t=T$.**

![Arbitrary COP system. Source: Author](figs/copDrawio.png){#fig:copSystem width=70% style="scale:1;"}
