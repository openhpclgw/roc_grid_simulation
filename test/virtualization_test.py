#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import (aggregate_current_vectors,
                            print_current_table,
                            energy_flow,
                            plot_heatmap,
                            plot_virtual_heatmap)


# I am using python 3.6.1

# assume square grid
# N = 100
# source = (10, 10, 5, 5)
# sink = (55, 10, 5, 50)
N = 40
source = (4, 4, 7, 10)
sink = [(20, 4, 1, 10),
        (10, 30, 2, 2)]
# source = (4, 4, 4, 12)
# sink = (12, 4, 4, 12)
cond_exp = -3
conductance = 10**cond_exp
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'
hp = HeatProblem(N, source, sink, conductance, src_val=10.)


m = ROCModel(mesh_size)

m.load_problem(hp)
m.run_spice_solver(virtualize=True)

plot_heatmap(m)
