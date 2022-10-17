# Literature Review

## Crashworthiness models: Masses-Springs-Dampers Systems
Crashworthiness models are models of vehicles used to analyze the safety of
their occupants in a crash. @Fig:crashworthiness shows an example. They can be
very elaborate, but very simple ones are also used on initial phases of vehicle
design. These simple ones are Masses-Springs-Dampers Systems (MSDS's), i.e.
they're comprised of only ideal masses, springs and dampers. They can be used to
estimate the maximum acceleration a driver would suffer in a crash.
Alternatively, we can look for values of the springs and dampers that would
minimize the maximum acceleration the drive would suffer in a crash.

![Example of Crashworthiness model. @Tbl:crashworthiness has the legend. Source: [@Mostafa2011-kc]](figs/crashworthiness.png){#fig:crashworthiness width=80%}


|  Mass No. |                    Vehicle components |
|:---------:|--------------------------------------:|
|          1|                    Engine and Radiator|
|          2|             Suspension and Front Rails|
|          3|             Engine Cradle and Shotguns|
|          4| Fire Wall and Part of Body on Its Back|
|          5|                               Occupant|
: Legend of @Fig:crashworthiness. Source: [@Mostafa2011-kc] {#tbl:crashworthiness}


TODO:

- Explain that all problems we optimize will be like crashworthiness optimization

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


