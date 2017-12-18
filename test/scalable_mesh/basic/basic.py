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


N = 40
source = (5, 5, 1, 1)
sink = [(35, 35, 1, 1)]
cond_exp = -3
conductance = 10**cond_exp
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'
hp = HeatProblem(N, source, sink, conductance, src_val=10.)

surface3d = False  # if True 3d surface plot is generated

m = ROCModel(mesh_size)

m.load_problem(hp)
m.run_spice_solver()
plot_heatmap(m, filename= 'basic_heat_meshsize{}'.format(mesh_size))
