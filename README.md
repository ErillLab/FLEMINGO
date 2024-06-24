![FLEMINGO Logo](https://github.com/ErillLab/FLEMINGO/images/FLEMINGO_logo.jpg)
## bacterial promoter modeling with chain-mode genetic programming

The idea of this project is to "learn" regulated promoter models, based on a set of user-provided sequences (referred to as the positive set), using a restricted Genetic Programming framework. The _restriction_ in this setup is that the models to be evolved are not open tree topologies (see [previous iteration](https://github.com/ErillLab/TF_GA) of the project), but chains connecting _recognizers_ via binary _connectors_.

The GP framework is tasked with evolving organisms that "bind" the positive dataset preferentially over a control dataset (referred to as the negative dataset). Optimal binding energies for the organisms on each sequence are computed using a variant of the Needleman–Wunsch algorithm.
