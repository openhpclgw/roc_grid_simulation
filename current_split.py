#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
from roc_model import ROCModel
import matplotlib.pyplot as plt
from heat_problem import HeatProblem
from analysis_utils import (run_current_split_analysis,
                            print_current_splits)

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
print_current_splits(run_current_split_analysis(m))
