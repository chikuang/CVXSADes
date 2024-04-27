# exact-design-SA

*[Chi-Kuang Yeh](https://chikuang.github.io/) & [Julie Zhou](https://www.uvic.ca/science/math-statistics/people/home/faculty/zhou_julie.php)*

*April 27, 2024*

---

### Description

Implementation of finding the optimal exact regression designs with complex settings using CVX in conjuction with the simulated annealing (SA) algorithm. The applications including 1) High-dimensional design, 2) Maximin design with competing optimal criteria and 3) Integer design problems.

### TODO
- [ ] It would be even better if we can use the symmetry property in the annealing algorithm for the cubic model. We can just work on the design points on [0, 1] and make it symmetric on [-1, 0].
- [ ] Cubic polynomial model; A- and D-optimality criteria
- [ ] Second-order linear models with 2 or 3 variables;  A- and I-optimality?
- [ ] Nonlinear models - 2-compartment model, peleg model, Emax model?
- [ ] GLMs with several design variables?
- [ ] A-optimality for exact design.

### Complete 
- In the D-optimal design criterion where the loss function is $-\log(\det(I(\xi)))$.
  - [x] Using an approximate design from CVX and N=21 for the discrete design space $S_N$, construct $n=20$ exact D-optimal design points in [-1, 1] via an annealing algorithm
  - [x] Using an approximate design from CVX and N=21 for the discrete design space $S_N$, construct $n=18$ exact D-optimal design points in [-1, 1] via an annealing algorithm.
  - [x] Using an approximate design from CVX and N=51 for the discrete design space $S_N$, construct $n=20$ exact D-optimal design points in [-1, 1] via an annealing algorithm.
  - [x] Using an approximate design from CVX and N=51 for the discrete design space $S_N$, construct $n=18$ exact D-optimal design points in [-1, 1] via an annealing algorithm.
