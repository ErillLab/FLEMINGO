![FLEMINGO Logo](https://github.com/ErillLab/FLEMINGO/blob/main/images/FLEMINGO_logo.jpg)

The idea of this project is to "learn" regulated promoter models, based on a set of user-provided sequences (referred to as the positive set), using an evolutionary algorithm. The _restriction_ in this setup is that the models to be evolved are not open tree topologies (see [previous iteration](https://github.com/ErillLab/TF_GA) of the project), but chains connecting _recognizers_ via binary _connectors_.

The algorithm is tasked with evolving organisms that "bind" the positive dataset preferentially over a control dataset (referred to as the negative dataset). Optimal binding energies for the organisms on each sequence are computed using a specific dynamic programming algorithm.
