# Example of max flow min cut for gender labeling problem

import maxflow
import numpy as np

#%%
d = np.array([179, 174, 182, 162, 175, 165]) # heights (data)
mu = [181, 165] # means of two classes
beta = 100 # weight of the prior term
w_s = (d-mu[0])**2 # source weight
w_t = (d-mu[1])**2 # sink weights
N = len(d) # number of graph nodes
indices = range(0,N) # an index for each person


# Create a graph with integer capacities.
g = maxflow.Graph[int](N,N)
# Add (non-terminal) nodes. Get the index .
nodes = g.add_nodes(N)
# Create edges between nodes
for i in range(0,N-1):
    g.add_edge(nodes[i], nodes[i+1], beta, beta)
# Set the capacities of the terminal edges.
for i in range(0,N):
    g.add_tedge(nodes[i], (d[i]-mu[1])**2, (d[i]-mu[0])**2)
# Run the max flow algorithm
flow = g.maxflow()
print('Maximum flow:', flow)

# displaying results
gend = 'MF'
for i in range(0,N):
    print('Person %d is estimated as %s' % (i, gend[g.get_segment(nodes[i])]))







