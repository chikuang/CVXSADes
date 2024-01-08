# exact-design-SA
Construct an exact design using CVX in conjunction with simulated annealing (SA).

## Authors:
1. [Julie Zhou](https://www.uvic.ca/science/math-statistics/people/home/faculty/zhou_julie.php)
2. [Chi-Kuang Yeh](https://chikuang.github.io/)


### TODO
- In the D-optimal design criterion where the loss function is $-\log(\det(I(\xi)))$.
  - [ ] Using an approximate design from CVX and N=21 for the discrete design space $S_N$, construct $n=20$ exact D-optimal design points in [-1, 1] via an annealing algorithm
  - [ ] Using an approximate design from CVX and N=21 for the discrete design space $S_N$, construct $n=18$ exact D-optimal design points in [-1, 1] via an annealing algorithm.
  - [ ] Using an approximate design from CVX and N=51 for the discrete design space $S_N$, construct $n=20$ exact D-optimal design points in [-1, 1] via an annealing algorithm.
  - [ ] Using an approximate design from CVX and N=51 for the discrete design space $S_N$, construct $n=18$ exact D-optimal design points in [-1, 1] via an annealing algorithm.
- Question: Can we use the symmetry property of the D-optimal design?
- Try A-optimality for the exact design afterward.
