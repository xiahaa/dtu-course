## simple FNN on 2D point set
![two layer fully connected neural network](../report/ex4/figures/nn1.png)

![initialization](../report/ex4/figures/ex2_1.png)
![initialization](../report/ex4/figures/ex2_2.png)
![initialization](../report/ex4/figures/ex2_4.png)

![initialization](../report/ex4/figures/ex1_1.png)
![initialization](../report/ex4/figures/ex1_2.png)
![initialization](../report/ex4/figures/ex1_4.png)

![initialization](../report/ex4/figures/ex3_1.png)
![initialization](../report/ex4/figures/ex3_2.png)
![initialization](../report/ex4/figures/ex3_4.png)

## CNN 
![KNN+VGG-16s features](../report/ex4/figures/ex5_1.png)
![retrained VGG-16](../report/ex4/figures/ex5_3.png)


## FNN on MNIST

1. 初始化有影响；
2. 加大加深网络不一定有作用，反而使得training变得很慢；
3. Adam，速度快一些，但是效果没有提升；
4. SGB+Momentum效果较好，效果取决于learning rate；
5. data normalization无用；
6. weight decay：无用；
7. data augmentation：有显著提升，0.5%；
8. minibatches: 获取的gradient更稳定；
9. early stopping: 有用，防止overfitting;

## Result:
Accuracy: 98.85, winner of 2019 MNIST competition for course Advanced Image Analysis
![misclassified images](https://github.com/xiahaa/cn/blob/master/report/ex4/figures/mnist1.png)

## TODO

1. CNN;

## git repo:

1. zzutk;
2. sunshine;
3. https://github.com/jakejhansen/Advanced_Image_Analysis_02503
