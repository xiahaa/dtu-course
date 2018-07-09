# Convex_Optimization
This repo is based on assignments done for course Convex Opimization.
It contains first-order method for large scale convex optimization problem.
It also contains my project assignment for course Convex Optimization where I did a semi-definite trajectory generation algorithm.

## add RRT*
Most RRT* algorithms you will find in github by searching RRT star contain several drawbacks:
1. almost none of them consider rewire which needs to iteratively refine the cost and its children's cost.[This one](https://github.com/olzhas/rrt_toolbox) rewires the costs of nearest neighbor but forgot to update their children's cost. [This one](https://github.com/saihv/rrtstar/blob/master/3D/RRTStar_3D.m) doesn't do any rewire. Basically, they can not be named RRT star, at most some kind of RRT variant.
2. I do not want to sample in a infinite space. Most time we will deal with would be a grid map, which is finite, so why should we sample infinitely???? In order to increase efficiency, I directly sample on the grid map and then mapping back to real metrics.
3. try to avoid loop as much as possible.

## add perlin noise
Map generator using perlin noise.

## lessons learned
1. free code doesn't mean free lunch!!!! Use you own judgement.
2. KD tree offered by VLFeat even slower than vectorized bruth-force search in MatLab, which seems to be strange...... I intended to use KD tree for nn search. However, after benchmark test, I found vectorized brute-force was even faster....So I abandoned to use KD tree..
