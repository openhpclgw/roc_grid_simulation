#! /usr/bin/env python
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import print_current_table

exp_prob_size = 5
prob_size = 2**exp_prob_size
source = (0, 0, 1, 1)
sink = (prob_size-1, prob_size-1, 1, 1)
exp_cond = -3
conductance = 10**exp_cond

filename='tmp/diag_v_{}'

mesh = ROCModel(prob_size)

def get_results(mesh, problem):
    mesh.load_problem(problem)
    mesh.run_spice_solver(filename.format(mesh.w))
    return mesh.final_grid

for v in range(1,5):
    print('Source potential = {} V'.format(v))
    hp = HeatProblem(prob_size, source, sink, conductance, src_val=v)
    res_grid = get_results(mesh, hp)
    print_current_table(mesh)
    print()
