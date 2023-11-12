# Methods {#sec:methods}

## Overview

The method used to achieve the objectives listed at @sec:objectives was to
implement a library that is capable of:

1. Defining arbitrary MSDSs, i.e. systems with arbitrary masses/springs/dampers
and also arbitrary initial/boundary conditions

2. Finding values of springs/dampers that minimize the maximum acceleration
that an arbitrary mass of an MSDS experiences after impact. The library must
implement two "flavors" of GA for this optimization: P-GA and E-GA.
- E-GA evaluates candidate solutions by explicitly integrating the MSDS in time.
- P-GA evaluates candidate solutions using PIMs which describe the position of
each mass as a function of *time* and of the values of the masses and springs
of the system

Using this library, we performed a few case studies in which we optimized MSDSs
using both P-GA and E-GA. The performance and the result of each algorithm were
then compared.