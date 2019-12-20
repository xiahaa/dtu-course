# Deep_Learning_in_Computer_Vision
Special course at DTU summer 19 - deep learning in computer vision

## Pre-reading
~~numpy: https://morvanzhou.github.io/tutorials/data-manipulation/np-pd/~~

~~matplotlib: https://matplotlib.org/tutorials/index.html~~

~~PIL or PILLOW: https://pillow.readthedocs.io/en/stable/handbook/tutorial.html~~ 貌似已经停止更新了，只支持很简单的一些图像操作，IO，基本几何变换，滤波。

torch

https://zhuanlan.zhihu.com/p/25572330

zhihu

https://www.zhihu.com/question/55720139

tutorial

https://pytorch.apachecn.org/docs/1.0/#/pytorch_with_examples

https://web.cs.ucdavis.edu/~yjlee/teaching/ecs289g-winter2018/Pytorch_Tutorial.pdf

https://pytorch.org/tutorials/

torch visualization

https://github.com/utkuozbulak/pytorch-cnn-visualizations/tree/master/src

https://discuss.pytorch.org/t/visualize-feature-map/29597/13

https://towardsdatascience.com/how-to-visualize-convolutional-features-in-40-lines-of-code-70b7d87b0030


## useful toolbox for optimization
1. Bayesian optimization: no need for gradient, hessian, jacobian. but treat f as a blackbox. 
GPyOpt: The tool for Bayesian Optimization: [example](https://github.com/ibalmeida/02901_Report) of using GPyOpt for automatically tune the neural network hyperparameters.

2. bayesian estimation needs to compute use the prior and likelihood to compute the maximum posterior estimation. if the likelihood is unknown or computationally intractable. can we approximate it? Via Learning Summary Statistic for Approximate Bayesian Computation via Deep Neural Network, an easy [example](https://github.com/SamuelWiqvist/project_02901_adv_ml/blob/master/project.ipynb).


