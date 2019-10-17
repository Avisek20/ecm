# Entropy c-Means Clustering

Entropy c-Means identifies fuzzy clusters by simultaneously optimizing two contradictory objectives: objective $f_1$ identifies compact clusters, whereas objective $f_2$ identifies overlapped clusters.

```
$\min f_1 = \sum_{i=1}^{N} \sum_{j=1}^{c} \mu_{ij} ||x_i - v_j||^2$

$\max f_2 = - \sum_{i=1}^{N} \sum_{j=1}^{c} \mu_{ij} \log(\mu_{ij})$
```

Paper Source: [Gupta A., Datta S. and Das S., "Fuzzy Clustering to Identify Clusters at Different Levels of Fuzziness: An Evolutionary Multiobjective Optimization Approach," in IEEE Transactions on Cybernetics, 2019.](https://ieeexplore.ieee.org/document/8692725)

## Implementations of Entropy c-Means

MATLAB implementations available in the **matlab** folder

Python implementations available in the **python** folder


