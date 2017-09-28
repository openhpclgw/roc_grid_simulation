#! /usr/bin/env python
import numpy as np
from codegen import *

# I am using python 3.6.1

# assume square grid
N = 4

# initialize grid
grid = np.zeros((N,N))

# heat up an arbitrary point
grid[2][2] = 1.

# heat up in gradient from the right
# for (i,j) in it.product(range(N), range(N)):
    # grid[i][j] = j*0.5

print(grid)
# initialize conductance table
# conductance = [[1. for i in range(N-1)] for j in range(N-1)]
conductance = 0.1

# print(grid)
m = ROCModel(4)
m.to_spice(grid, conductance)
