#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
from roc_model import ROCModel
import matplotlib.pyplot as plt
from heat_problem import HeatProblem
from analysis_utils import (aggregate_current_vectors,
                            print_current_table,
                            energy_flow,
                            node_potentials,
                            nodal_current_tup)


# I am using python 3.6.1

# assume square grid
# N = 100
# source = (10, 10, 5, 5)
# sink = (55, 10, 5, 50)
N = 20
source = (0, 0, 1, 1)
sink = (19, 19, 1, 1)
# source = (4, 4, 4, 12)
# sink = (12, 4, 4, 12)
cond_exp = -3
conductance = 10**cond_exp
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'
hp = HeatProblem(N, source, sink, conductance, src_val=10.)


m = ROCModel(mesh_size)

m.load_problem(hp)
m.run_spice_solver(hp)


# record initial current splits
init_nodal_current_splits = [[nodal_current_tup(m.nodes[i][j])
        for i in range(mesh_size)] for j in range(mesh_size)]

final_nodal_current_splits = [[{}
        for i in range(mesh_size)] for j in range(mesh_size)]

# print(init_nodal_current_splits)

def split_dif(biased, base):
    return {d:abs(biased[d]-base[d]) for d in ('E', 'W', 'N', 'S')}

def mask_split(split, mask_key):
    split[mask_key] = 0.0

def normalize_split(split):
    val_sum = sum([v for _,v in split.items()])
    if val_sum != 0:
        return {k:v/val_sum for k,v in split.items()}
    else:
        return split

# apply bias and measure the deltas
for i,j in it.product(range(mesh_size), range(mesh_size)):
    node = m.nodes[i][j]
    base_split = init_nodal_current_splits[i][j]
    fin_splits = final_nodal_current_splits[i][j]

    for d,a in node.ammeters.items():
        sys.stdout.write('\rRunning bias: {}, {}'.format(str((i,j)), d))
        a.set_bias(1)
        m.run_spice_solver(hp)
        biased_split = nodal_current_tup(node)
        a.reset_bias()
        dif_split = split_dif(biased_split, base_split)
        mask_split(dif_split, d)
        fin_split = normalize_split(dif_split)
        fin_splits[d] = fin_split

# print the final table
row_format = '{}\t\t{}\n\t\t{}\n\t\t{}\n\t\t{}'
dict_format = 'E:{E:0<.2f} W:{W:0<.2f} N:{N:0<.2f} S:{S:0<.2f}'
for i,j in it.product(range(mesh_size), range(mesh_size)):
    tmp_splits = final_nodal_current_splits[i][j]
    # print(*tmp_splits)
    print(row_format.format(str((i,j)),
        *[dict_format.format(**d) for _,d in tmp_splits.items()]))

    
    print()

