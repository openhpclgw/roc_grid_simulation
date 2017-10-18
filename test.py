#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
from roc_model import ROCModel
import matplotlib.pyplot as plt

class HeatProblem(object):
    # bbox : (left, top, width, height)
    def __init__(self, N, source, sink, conductance, src_val=1., sink_val=0.):
        self.N = N
        self.source = source
        self.sink = sink
        self.conductance = conductance
        self.src_val = 1.
        self.sink_val = 0.

    def __iter_bbox(self, bbox):
        for i,j in it.product(range(bbox[1], bbox[1]+bbox[3]),
                              range(bbox[0], bbox[0]+bbox[2])):
            yield (i,j)

    def source_iter(self):
        for idx in self.__iter_bbox(self.source):
            yield idx

    def sink_iter(self):
        for idx in self.__iter_bbox(self.source):
            yield idx

    def __is_in_bbox(self, bbox, point):
        i, j = point
        left, top, width, height = bbox

        if j < left or j >= left+width:
            return False
        if i < top or i >= top+height:
            return False

        return True

    def is_source(self, p):
        return self.__is_in_bbox(self.source, p)

    def is_sink(self, p):
        return self.__is_in_bbox(self.sink, p)

    def gen_matrix(self):
        mat = np.zeros((N,N))
        print(self.source)
        for (i,j) in self.__iter_bbox(self.source):
            mat[i][j] = self.src_val
        for (i,j) in self.__iter_bbox(self.sink):
            mat[i][j] = self.sink_val

        return mat
        


# I am using python 3.6.1

# assume square grid
N = 100
source = (0, 70, 5, 29)
sink = (70, 20, 29, 5)
cond_exp = -3
conductance = 10**cond_exp
num_iters = 1
pos = 1
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'
hp = HeatProblem(N, source, sink, conductance, src_val=100.)


# sum_result = np.zeros((mesh_size,mesh_size))
m = ROCModel(mesh_size)

# for i in range(num_iters):
    # result = m.run_spice_solver()
    # for s,r in zip(sum_result, result):
        # s += r

debug = False
if debug:
    for i in range(N):
        for j in range(N):
            print('{0: >11.5f} '.format((result[i][j])*10.**15), end='')
        print()

# save the input grid as well
# plt.imshow(grid, cmap='hot', interpolation='nearest')
# plt.savefig('tmp')

def cross_plot(i, j):
    data_x = [sum_result[i][idx] for idx in range(mesh_size)]
    data_y = [sum_result[idx][j] for idx in range(mesh_size)]
    data_cross = [sum_result[idx][idx] for idx in range(mesh_size)]
    figs, axes = plt.subplots(3,3)
    axes[2][1].plot(range(mesh_size), data_x)
    axes[1][2].plot(range(mesh_size), data_y)
    axes[0][0].plot(range(mesh_size), data_cross)
    axes[1][1].imshow(sum_result, cmap='hot', interpolation='nearest')

# def heat_plot(data, hp):
    # zero out the source and the sink
    # for i,j in hp.source_iter():
        # data[i][j] = 0.
    # for i,j in hp.sink_iter():
        # data[i][j] = 0.

    # plt.imshow(data, cmap='hot', interpolation='nearest')

def numerical_solve(hp, num_steps):
    N = hp.N
    # grid = hp.gen_matrix()
    # grid2 = hp.gen_matrix()
    grids = (hp.gen_matrix(), hp.gen_matrix())
    c = hp.conductance

    for step in range(num_steps):
        for (i,j) in it.product(range(N), range(N)):
            if hp.is_source((i,j)):
                continue
            if hp.is_sink((i,j)):
                continue
            ing=step%2
            outg=1-ing
            tmp_sum = 0.
            if i-1 >= 0:
                tmp_sum += grids[ing][i-1][j]
            if i+1 < N:
                tmp_sum += grids[ing][i+1][j]
            if j-1 >= 0:
                tmp_sum += grids[ing][i][j-1]
            if j+1 < N:
                tmp_sum += grids[ing][i][j+1]
            grids[outg][i][j] = (tmp_sum)/4.

        abs_delta = 0.
        for (i,j) in it.product(range(N), range(N)):
            abs_delta += abs(grids[0][i][j]-grids[1][i][j])
        print(abs_delta)
        if abs_delta < 0.01:
            break

    return grids[num_steps%2]

# def quiver_plot(U,V):

# tmp_x_loc = int(pos/int(N/mesh_size))
# tmp_y_loc = int(pos/int(N/mesh_size))
# cross_plot(tmp_x_loc,tmp_y_loc)
# plt.savefig(img_name.format(gr_sz=N, ms_sz=mesh_size))
# heat_plot(numerical_solve(hp,10000), hp)
results, U, V = m.run_spice_solver(hp)

fig, axes = plt.subplots(1,1)
axes.imshow(results, cmap='hot', interpolation='nearest')
# print([i for i in range(mesh_size)])
# print([i for i in range(mesh_size-1, -1, -1)])

axes.streamplot(np.array([i for i in range(mesh_size)]),
               np.array([i for i in range(mesh_size)]),
               U, V, color='green', linewidth=1, density=0.5)

# heat_plot(results, hp)
plt.show()
