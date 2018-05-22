#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import (aggregate_current_vectors,
                            print_current_table,
                            energy_flow,
                            plot_surface,
                            plot_heatmap)


# I am using python 3.6.1

# assume square grid
# N = 100
# source = (10, 10, 5, 5)
# sink = (55, 10, 5, 50)

#Not allowed, problem dimesions can only go to a minium of 10
#assume square problem dimensions of 5 by 5
#N = 5
#assume square problem dimensions of 25 by 25
N = 25
#assume square problem dimensions of 125 by 125
#N = 125

#Grid Layout
#y

#0
#1
#2
#3
#.
#.
#.
#24
#   0   1   2   ...  24 x


# Is the tuple laid onto the problem, or the mesh(grid) => The Problem, Good!
# Setup for the 25 problem, does not allow for variable source or sink => Bad
# Where is the source value set?
# I belive the source and sinks can not overlap => WRONG, Sources and Sinks CAN overlap
# I belive the source and sinks can not go out of bounds of the problem=> yes
# Index starts at 0 => Yes

# I am not going to overlap becuase I dont know how it effects values? =>
# The tuple is in form (left_pos(x axis), top_pos(y axis), width(x axis), height(y axis))
#top
source = (0, 0, 25, 1)
#left, bottom, right
sink = [(0, 1, 1, 24),(0, 24, 24, 1),(24, 1, 1, 24)]

# sink = [(11, 2, 1, 10), (5, 15, 2, 2)]
# source = (4, 4, 4, 12)

cond_exp = -3
conductance = 10**cond_exp

#mesh dimension entered by user
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'
hp = HeatProblem(N, source, sink, conductance, src_val=10.0)

surface3d = False  # if True 3d surface plot is generated

m = ROCModel(mesh_size)

m.load_problem(hp)
m.run_spice_solver()

eflow_data = energy_flow(m)
print('Source in : {}'.format(eflow_data['src_in']))
print('Source out: {}'.format(eflow_data['src_out']))
print('Sink in   : {}'.format(eflow_data['snk_in']))
print('Sink out  : {}'.format(eflow_data['snk_out']))

print_current_table(m)

# current flows currently doesn't work correctly with the new grid
# logic. It would probably require some post processing on the grid
# itself (which sucks)
if not surface3d:
    plot_heatmap(m)
else:
    plot_surface(m)
