#! /usr/bin/env python
import sys
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import *

size = int(sys.argv[1])
N = size
mesh_size = int(sys.argv[2])
source = (0, 0, 1, size)
sink = [(1, 0, size-1, 1), (size-1, 0, 1, size-1), (1, size-1, size-1, 1)]
conductance = 10**-3  # this'll be used as resistance directly
hp = HeatProblem(N, source, sink, conductance, src_val=286.)

m = ROCModel(mesh_size)
m.load_problem(hp)
m.run_spice_solver(filename='problem')

from_mesh = m.final_grid
from_comsol = load_grid_from_comsol_csv('../comsol_data/test.csv')
plot_errmap(from_mesh, from_comsol)
