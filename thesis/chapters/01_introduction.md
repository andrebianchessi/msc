# Introduction
TODO: write introduction with the following key points

- Explain relevance of mechanical optimizations.
- Show literature examples of metamodels used for genetic algorithm optimization.
- Explain that metamodels require a high cost for their training, but might
accelerate the overall optimization for, once trained, the AI is faster than
the FEM model. That, however, is not trivial/straightforward. It might "not be worth it".
- Explain that data fit models are only capable of interpolations, and not extrapolations
- Introduce PINNs, citing that they've been used in studies in many areas (CFD, topology opt etc).
- Explain that PINNs are new and promising models because:
    - Incorporate physics knowledge into the model, which might increase the extrapolation capability
    - Can be trained without labeled data; thus don't require running many FEM simulations beforehand.
- Explain that this study focuses on exploring the usage of PINNs for genetic algorithm optimization.
We set out to analyze how much performance gain they can provide when compared to traditional methods (FEM, discrete simulations),
how that changes with respect to the system's complexity, and how they perform if used to simulate conditions outside of the ones they were trained on,
which could allow models to be reused on subsequent optimizations of the same system with different initial conditions.
- Highlight that our focus is on exploring Genetic Algorithm + PINN optimization, and not necessarily on assessing the benefits of applying this technique to a specific problem. We know that results to some case studies are not trivially generalizable, but they may still increase our understanding on potentials and limitations of this promising new technique, provide some insights and new ideas worth pursuing, and facilitate future work that might be similar to this one.
- Inspired by lumped crashworthiness models, we chose to optimize mass/spring/damper systems. They can easily be made arbitrarily simple/complex and require a small implementation overhead, which allows us to spend less effort on implementing FEM solvers for example.