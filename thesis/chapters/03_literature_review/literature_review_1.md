## Discrete element method {#sec:dem}

In order to explicitly integrate a mechanical system with respect to time, the
**ODE that describes the system** is necessary. The discrete element method is
a method to obtain the ODE for [CMs](#sec:cms).

### Local Matrices
@Fig:discreteElementSimple shows the simplest system that contains ideal masses,
springs and dampers. The spring has its natural/relaxed length when $x_0=x_1=0$.

![CM of 2 masses, 1 spring and 1 damper. Source: Author](figs/discreteElementSimple.png){#fig:discreteElementSimple width=80% style="scale:1;"}

Considering that the elastic force $F_k$ is linear with respect to the displacement ($F_k = k \cdot x$) and that
the damping force $F_c$ is linear with respect to the speed ($F_x = c \cdot \dot{x}$),
from Newton's second law we obtain @eq:demLocalEq, in which $\dot{x_i}$ and $\ddot{x_i}$
represent the first and second time derivative
of the displacement (with respect to the springs' relaxed position) of the mass $i$.

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

### Assembling the global matrix

By adding another mass, spring and damper to @fig:discreteElementSimple, we can obtain the system at @fig:discreteElementSimple2.

![CM of 2 masses, 2 spring and 2 dampers. Source: Author](figs/discreteElementSimple2.png){#fig:discreteElementSimple2 width=80% style="scale:1;"}

Equations of motion of system illustrated at @fig:discreteElementSimple2, also obtained by Newton's second law, are:

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

![CM of 3 masses, 3 spring and 3 dampers. Source: Author](figs/discreteElementSimple3.png){#fig:discreteElementSimple3 width=80% style="scale:1;"}

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

![CM of 3 masses, 4 spring and 4 dampers. Source: Author](figs/discreteElementSimple4.png){#fig:discreteElementSimple4 width=80% style="scale:1;"}

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
global matrices from the elements' local matrices. The reader can find in more
detail at [@Logan2007-bq; @Alves2020-gz].
