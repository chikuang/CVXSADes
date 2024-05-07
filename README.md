# CVXSADes: CVX-simulated annealing algorithm for finding exact regression designs

*[Chi-Kuang Yeh](https://chikuang.github.io/) & [Julie Zhou](https://www.uvic.ca/science/math-statistics/people/home/faculty/zhou_julie.php)*

*May 05, 2024*

---

### Description

Implementation of finding the optimal exact regression designs with complex settings using CVX in conjunction with the simulated annealing (SA) algorithm. The applications include 1) high-dimensional design, 2) Maximin design with competing optimal criteria, and 3) integer design problems.

### Computing Environment

Windows Machine with an i9 processor equipped with MATLAB 2023 version B update 2 and CVX version 2.2.

### Applications and models

1. Two-variables logistic regression model with interaction.
    - This model was used in [Haines et al.](https://linkinghub.elsevier.com/retrieve/pii/S0378375817301441) (2018) to study the usage of two insecticides, rotenone and pyrethrin. They were sprayed either alone or jointly to a species of fly, the chrysanthemum aphid. The response of interest was the number of insects dead, and the explanatory variables were the dose concentrations of the two insecticides.
2. Seven-dimensional logistic regression model with some interactions. 
   - This model was used in [Xu et al.](https://ieeexplore.ieee.org/document/8598838/) (2019). This paper concerns a car refuelling experiment with discrete and continuous factors and selected pairwise interactions.
3. Group design for disease prevalence.
   - This model was used in [Huang et al.](https://academic.oup.com/jrsssb/article/79/5/1547/7041031?login=false) (2017). They study the prevalence of a train using a test with uncertain sensitivity and specificity. A group testing study often aims to estimate the prevalence of a rare disease or a particular trait. Group testing is frequently used in studies where testing individuals for a trait is expensive and individual samples are relatively plentiful. 
4. Maximin design with three competing design criteria.
   - This model was used in [Bretz et al.](https://onlinelibrary.wiley.com/doi/10.1002/sim.3802) (2010) and [Wong and Zhou](https://www.tandfonline.com/doi/full/10.1080/10618600.2022.2104858) (2023) for a dose-finding study for an anti-asthmatic drug. There are 4 competing models, namely 1) Linear, 2) EMax 1, 3) Emax 2 and 4) Logistic regression models. When to use which model is unclear as each model possesses their own advantages, thus, this study aims to find the optimal design that maximizes the minimum efficiency of those four competing models.

### Reference paper
Yeh, C.-K. & Zhou, J. (2024+). CVXSADes: A stochastic algorithm for constructing optimal exact regression designs with single or multiple objectives. [arXiv preprint](https://arxiv.org/abs/2405.02983), [ReseachGate link](https://www.researchgate.net/publication/380295508_CVXSADes_a_stochastic_algorithm_for_constructing_optimal_exact_regression_designs_with_single_or_multiple_objectives)
