#! /usr/bin/env python
import sys
from roc_model import ROCModel
import matplotlib.pyplot as plt
from heat_problem import HeatProblem
from analysis_utils import *
from att_sweep_helper import *

def resistance(size, n1, n2):
    # second 0.25 percentile from top
    if ((n1[0] < size/2 and n1[0] >= size/4) or
        (n2[0] < size/2 and n2[0] >= size/4)):
        # first 0.25 percentile from left
        if n1[1] < size/4 or n2[1] < size/4:
            return 0.1*multiplier
        # last 0.25 percentile from right
        if n1[1] >= size-(size/4) or n2[1] >= size-(size/4):
            return 100000*multiplier
    return 0.001*multiplier


size = mesh_size
source = (0, 0, 1, int(size/4))
sink = (size-1, size-int(size/4), 1, int(size/4))
hp = HeatProblem(size, source, sink, resistance)

working_dir = 'test/attenuation_sweep/case3/size'+str(mesh_size)+'/'
out_filename = working_dir+'att'+str(base_att)+'/{}_'+str(size)
scr_name = working_dir+'opt_case3_'+str(size)+'_att'+str(base_att)

att_sweep(hp, scr_name, out_filename, working_dir)
