## Physics-Informed Machine Learning {#sec:pim}

We start by explaining the basic concepts of Machine Learning, and then show how
the same ideas, with some small modifications, are applied to solve problems
which are described by partial differential equations.

The problems we choose to solve are very simple, but they're great for illustrating
the basic principles.

### Introduction to Machine Learning

Consider the data in @tbl:mlExample.

|  $t$ | $y$ |
|:----:|:---:|
|    10|  1.0|
|    20|  2.1|
|    30|  3.0|
: Example dataset {#tbl:mlExample}

Let's say we want to model the relationship between $x$ and $y$ in a linear
model in the form $\bar{y}(t) = at + b$, where $a$ and $b$ are real constants.
This would allow us to approximate the values of $y$ for values of $x$ other
than the ones we have in our dataset. Note that the model we choose is
arbitrary. We could pick any other much more complicated model of higher order,
but we'll stick to this linear one because the principles are the same.

For example, a possible model would be $\bar{y}(t) = 2t + 5$. That would, however, be terrible to model our data, as the difference between the model and the real values are very big. Hence, we want to find values of $a$ and $b$ that cause our model to make predictions close to the real values. Ideally, we'd want $\bar{y}(t_i) = y_i$ for all $t_i$ from @tbl:mlExample.
An expression that we can use to measure how well our model fits to the data is to calculate the sum of the square of the residues:

$$
R = \sum_{t_i} (\bar{y}(t_i)-y_i)^2 = (\bar{y}(10)-y_0)^2 + (\bar{y}(20)-y_1)^2 + (\bar{y}(30)-y_2)^2
$$

We take the square of the residues so that residues of opposite signs do not
cancel each other out. This is equivalent to taking the absolute value of the
residue, but much more mathematically convenient.

Replacing the definition of $\bar{y}$ and the $y_i$ values from @tbl:mlExample:

$$
R = (10a + b -1)^2 + (20a + b -2.1)^2 + (30a + b - 3.0)^2
$$

From now on, our problem becomes finding values of $a$ and $b$ that minimize $R$. 

We can find those values with the gradient descent method.
Since we can calculate the derivatives of $R$ with respect to $a$ and $b$, we can also calculate its gradient, which can be thought of as a vector in the $a,b$ plane which points in the direction in which $R$ increases the most:

$$
\nabla R =
\begin{bmatrix}
\frac{\partial R}{\partial a} \\
\frac{\partial R}{\partial b}
\end{bmatrix}
=
\begin{bmatrix}
2800 a + 120 b - 284 \\
120 a + 6 b - 12.2
\end{bmatrix}
$$

The basic idea in the gradient descent method is to start with random $a$ and $b$ values, and make successive steps in the opposite direction of the gradient. The code listing bellow, written in python, exemplifies how gradient descent can be performed:

```
def grad(a,b):
    return (2800*a + 120*b - 284, 120*a + 6*b - 12.2)

a = 1
b = 1
l = 0.00001
for i in range (10000000):
    g = grad(a,b)
    a = a - l*g[0]
    b = b - l*g[1]

print(a,b)
```

The listing above outputs $a = 0.0999$ and $b = 0.0333$, which is, up to numerical precision, exactly like the analytical solution.

Thus, we find that our model has the following expression, which yields results very close to @tbl:mlExample:

$$
\bar{y}(t) = 0.099t + 0.033 
$$

Since we're just using an example to illustrate the principles of ML, we chose a very simple model, with a single input, and very few data points. However, the general workflow is the same for real like applications:

1. Define a model
2. Choose loss function that has the model's parameters (in our example those were just $a$ and $b$) as arguments, and measures, using the data we have available, how well our model's predictions match the expected outputs
3. Minimize the loss function to find the optimal parameters for the model

### Introduction to Physics-Informed Machine Learning

What if instead of having values for $t$ and $y$, we had an expression for $\ddot{y}$, the second time derivative of $y$?
For example: Consider a ball that is thrown upwards with an initial velocity of $10$m/s. With $y$ as the height of the ball and considering $10$m/s$^2$ as gravity's acceleration, the physical equation that governs this motion is:

$$
\ddot{y} = -10
$$
{#eq:pimExample}

We can solve this problem in the same way as we did the previous one by just doing some modifications.
First, we must define the model we want to use to approximate $y$.
Let $\bar{y}(t)$ be the model, given by:

$$
\bar{y}(t) = a\cdot t^2 + b\cdot t + c
$$
{#eq:pimModel}

We chose this because we can then compare the solution we obtain with the analytical solution $y(t) = -5t^2+10t$.

The next step is to define a loss function. In this case, we only have value for $y$ at $t=0$ (considering that the ball starts at position $y=0$). We also have value for the first derivative of $y$ at $t=0$, which is the initial speed. For all other instants of time, we only know the second derivative of $y$ from @eq:pimExample.
To write a loss function, we can first define multiple values of $t$ which are of interest. Then, the loss function can calculate the residues for $y$ and $\dot{y}$ at $t=0$, and for $\ddot{y}$ for all other time instants:

Let $t = [0, 0.1, 0.2, ... , 1.0]$. The loss function $R$ is given by:

$$
R = (\bar{y}(0)-y|_{t=0})^2 + (\dot{\bar{y}}(0)-\dot{y}|_{t=0})^2 +
\sum_{i=0}^{10} (\ddot{\bar{y}}(0.1 i)-\ddot{y})^2
$$

Note that we can use this technique because we chose what we want the model $\bar{y}$ to be like; so we can easily calculate its time derivative. Replacing $\bar{y}(t)$ with its definition from @eq:pimModel, $y|_{t=0}$ with $0$ and $\dot{y}|_{t=0}$ with $10$ (the initial velocity), we get:

$$
R = (c-0)^2 + (b-10)^2 +
\sum_{i=0}^{10} (2a+10)^2
$$

All we're left to do now is to find the values of $a$, $b$ and $c$ that minimize $R$.
The gradient of R is given by:
$$
\nabla R =
\begin{bmatrix}
\sum_{i=0}^{10} 4(2a+10) \\
2(b-10) \\
2c
\end{bmatrix}
$$

A simple python script that can be used to perform gradient descent is the following:

```
nT = 10

def grad(a,b,c):
    ga = 0
    for i in range (nT+1):
        t = 1/nT*i
        ga += 4*(2*a+10)
    
    gb = 2*(b-10)

    gc = 2*c
    
    return (ga,gb,gc)

a = 0
b = 0
c = 0
l = 0.0001
for i in range (100000):
    g = grad(a,b,c)
    a = a - l*g[0]
    b = b - l*g[1]
    c = c - l*g[2]
    print(a,b,c)
```

The output of the above script is ```-4.9999 9.9999 0.0``` which is, up to
numerical precision, exactly like the analytical solution.

@Thuerey2021-ut is an excellent source that covers, in great depth and breadth, the fundamentals of physics based models and the usage of more sophisticated models such as deep neural networks.