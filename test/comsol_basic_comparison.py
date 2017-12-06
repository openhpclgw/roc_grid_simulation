#! /usr/bin/env python
import sys
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import *


# in COMSOL I still can't create a 2D region that is to be taken as
# fixed temperature in our case. Which may not make sense with
# physical/realistic use cases, anyway.
# So, here we are creating a heat problem that includes source and sink
# regions, but then when comparing against comsol data, we need to get
# the center slice which doesn't have any source or sink within.
# so we need to create a heat problem larger than the comsol data, so
# that we can get the center slice that has the same data resolution
size = int(sys.argv[1])
mesh_size = int(sys.argv[2])

comp_size = int(size/(mesh_size-2))
if comp_size != size/(mesh_size-2):
    print('Divisibility error')

prb_size = size+2*comp_size  # one component on each side

source = (0, 0, comp_size, prb_size)
sink = [(comp_size, 0, prb_size-comp_size, comp_size),
        (prb_size-comp_size, 0, comp_size, prb_size),
        (comp_size, prb_size-comp_size, prb_size-comp_size, comp_size)]


conductance = 10**-3  # this'll be used as resistance directly
hp = HeatProblem(prb_size, source, sink, conductance, src_val=1.)

m = ROCModel(mesh_size)
m.load_problem(hp)
m.run_spice_solver(filename='problem')

from_mesh = m.final_grid[comp_size:comp_size+size,
                         comp_size:comp_size+size]
from_comsol = load_grid_from_comsol_csv('comsol_data/left_src_others_sink_comsol.csv')

# print_error_table(from_mesh, from_comsol)
plot_errmap(from_mesh, from_comsol)


# print_node_potentials(m)
# plot_heatmap(m)
