#! /usr/bin/env python
from codegen import *

# I am using python 3.6.1

# assume square grid
N = 4

# initialize grid
grid = [[0. for i in range(N)] for j in range(N)]

# heat up an arbitrary point
grid[2][3] = 1.

# initialize conductance table
conductance = [[1. for i in range(N-1)] for j in range(N-1)]

print(grid)
m = SpiceModel(grid, conductance)
m.to_spice()
