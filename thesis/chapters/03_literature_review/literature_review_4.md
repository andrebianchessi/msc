## Meta-models for mechanical optimization

[@Gu2018-uk; @Gu2018-tf; @Lee2022-uz] are among the main sources of inspiration
for this work. They all studied the usage of metamodels in conjunction with
Genetic Algorithm for optimization of materials/solid structures in a vast
design space, and were very successful. The metamodels used provided a very
significant improvement in performance, and the solutions that were obtained
were very efficient. [@Driemeier_undated-za] studied using metamodels to predict
the behavior of highly non-linear 3D lattice structure.

All the studies mentioned above had to go through a very computationally
expensive process of generating labeled data and training the metamodels.

[@Zhang2021-vj; @Zhang2022-jp] have studied PIMs, which can be trained without labeled data. These models make it possible to *bypass* the stage of training data generation. This led us to wonder if PIMs might be the ideal metamodels for metamodel-based-GAs. However, we didn't find, in the literature, broad studies of this technique itself, i.e. studies in which problems with different complexities were analyzed.