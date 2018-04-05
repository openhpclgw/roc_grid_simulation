#! /usr/bin/env python
import sys
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import *

def resistance(size, n1, n2):
    return 0.001

size = int(sys.argv[1])

N = size
mesh_size = size
source = (0, 0, 1, size)
sink = [(size-1, 0, 1, size)]
hp = HeatProblem(N, source, sink, resistance)

gen_interconnect_script = False
# assert not gen_interconnect_script


m1 = ROCModel(mesh_size)
m1.load_problem(hp)
m1.run_interconnect_solver(filename='interconnect_'+str(size),
                          gen_script=gen_interconnect_script,
                          get_results=not gen_interconnect_script)

# plot_heatmap(m1, current_flow_plot=None)

if not gen_interconnect_script:
    m2 = ROCModel(mesh_size, sidelinks=False)
    m2.load_problem(hp)
    m2.run_spice_solver()

    normed_m1_grid = normalize_grid(m1.final_grid)
    normed_m2_grid = normalize_grid(m2.final_grid)

    # print(normalize_grid(m1.final_grid))
    # print(m2.final_grid)

    plot_errmap(normed_m1_grid, normed_m2_grid,
                filename='test/optical_comparison/opt_minus_elec')

    max_err = (normed_m1_grid-normed_m2_grid).max()
    print(normed_m1_grid)
    print(normed_m2_grid)
    print(max_err)
