# bayesian-scientifical-computing
This repository contains the lecture notes and exercises of Bayesian Scientifical Computing course taught by Daniela and Erkki at DTU 2019.

## bayesian scientifical computing

multivatiant:
$$
\mathbf{X} = \phi(\mathbf{Y})
$$

$$
\int_{\mathbf{B}} \pi_{\mathbf{X}} (x) dx=\int_{\phi^{-1}\mathbf{B}} \pi_{\mathbf{X}}(\phi (\mathbf{Y})) |J| dy
$$
where $J$ is the jacobian of $\phi$ over $y$.

**well posed**
1. solution exists
2. solution is unique
3. small perturbabtion on the variables results in small perturbation on the solution (stability)

ill posed is the opposite of well posed. if not well posed, then ill posed.

**span**
the usage of span is to find the minimal representation and eliminate redundant variables.

# summary of bayesian scientifical computing
## day1
page 70, $Ax=b$ exist of solution.
page 73, range, null space.
page 85, condition number.
page 91, null space means loss of information.
page 95, 	four fundamental spaces.

## day2: prior and likelihood
page 23: smoothness prior, first order and second order.
page 29: Whittle-Matern prior, $\lambda$ > 0 is the correlation length, $\beta$ is the smoothness order.
page 31: how to add struture.
page 41: **data driven prior.**
page 45-51: how to draw from probability distribution. actually draw from cumulative probability density. the problem is how you inverse from a probability to a variable.
page 55: Tikhonov regularized solution.

exercise: random draw
1D

![img](https://github.com/xiahaa/bayesian-scientifical-computing/blob/master/docs/fig1.png)

2D

![img](https://github.com/xiahaa/bayesian-scientifical-computing/blob/master/docs/fig3.png)

adding structure

![img](https://github.com/xiahaa/bayesian-scientifical-computing/blob/master/docs/fig4.png)

hypermodel

![img](https://github.com/xiahaa/bayesian-scientifical-computing/blob/master/docs/fig5.png)



## day3: conditional 
page 19-22: schur complement and its relationship to conditional gaussion.
page 28-29: apply schur complement to linear inverse problem and show the posterior is always more informative than the prior.
page 27, 33 are equivelent thanks to the Sherman-Morrison-Woodbury identity (aka matrix inverse lemma). proof is given 34-35.
exercise: **Kriging** interpolation, is actually a Maximum A Posterior estimation.

![img](https://github.com/xiahaa/bayesian-scientifical-computing/blob/master/docs/fig6.png)


## day 4: solver of Ax=b
page 6: conjugate gradient method;
page 8-9: generalization of cgls;
page 14: MAP via CGLS solver;
page 12: early stopping using Morozov Discrepancy Principle (MDP).
page 17: prior conditioner.
ex1: gaussian process simulation, note that there is a analytical solution.

![img](https://github.com/xiahaa/bayesian-scientifical-computing/blob/master/docs/fig7.png)

ex2: deconvolution using bayesian MAP estimation. how to use IAS solver.

![img](https://github.com/xiahaa/bayesian-scientifical-computing/blob/master/docs/fig8.png)


## day 5: dynamic inverse problem
page 3-4: definition of dynamic inverse problem.
page 36: **Sampling Importance Resampling (SIR)** based particle filtering.
page 37-40: use predictor and auxliary particle to improve normal particle filtering. 
**Survival of the Fittest (SOF)**:
1. 采样，assign权重；
2. propagate, 此时不加noise;
3. 算似然，归一化，重要性采样；
4. 在重采样之后的粒子上，加noise
5. 重新算扰动之后的似然，更新粒子的weight；

**Estimating parameters: Sequential Monte Carlo**
重点在42，其余的具体实施可比照SOF。

**Ensemble Kalman filter**
page 58-60:
1. 采样；
2. evolution+扰动；
3. 算mean cov;
4. 扰动观测；
5. 求解MAP更新粒子；
6. 算后验mean cov;

比粒子滤波需要的粒子数目大大减少，相比于EKF，对非线性的适应性可能更强。增加的计算量在于求解一个优化为题。
EnKF的名字的来源在于，如果观测模型是线性，则MAP等同于一个KF的update的过程，那么EnKF相当于对每一个粒子做一个KF，再综合。

page 62: EnKF + Add parameters.

![img](https://github.com/xiahaa/bayesian-scientifical-computing/blob/master/docs/SOS1.png)
![img](https://github.com/xiahaa/bayesian-scientifical-computing/blob/master/docs/SOS2.png)
![img](https://github.com/xiahaa/bayesian-scientifical-computing/blob/master/docs/SOS3.png)

