## Explicit Time Integration (ETI) {#sec:eti}

The [Discrete Element Method](#sec:dem) can be use to obtain a [CM](#sec:cms)'s
ODE, i.e. the system's state vector $X$ and an expression $\dot{X}(X)$ that
calculates the state vector's time derivative based on the current state. The
ODE can then be integrated with an explicit method such as Forward Euler.
Starting with the initial conditions, which are given, these methods recursively
calculate the state of a system at a later time based on its state at a current
time.

**The output of the integration is the time response of the system, i.e. the
displacement, velocity and acceleration of all the masses. This way, the maximum
acceleration that each mass will experience - which is the loss function that is
to be minimized with the [COPs](#sec:cop) - can be obtained**.

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
\begin{aligned}
&
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
\\\\
&
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
\end{aligned}
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

Thus, the ODE is defined as:
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

If the ODE were to be integrated with a simple Forward Euler using a $0.1$ time-step, the first time-step would be:

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

### Software

See @sec:software_eti for how we implemented explicit
time integration of [CMs](#sec:cms).
