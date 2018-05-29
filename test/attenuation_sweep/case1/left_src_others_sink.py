from heat_problem import HeatProblem
from att_sweep_helper import *

def resistance(size, n1, n2):
    return 0.001*multiplier

size = mesh_size
source = (0, 0, 1, size)
sink = [(1, 0, size-1, 1),
        (size-1, 0, 1, size-1),
        (1, size-1, size-1, 1)]
hp = HeatProblem(size, source, sink, resistance)

working_dir = 'test/attenuation_sweep/case1/size'+str(mesh_size)+'/'
scr_name = working_dir+'opt_case1_'+str(size)+'_att'+str(base_att)
out_filename = working_dir+'att'+str(base_att)+'/{}_'+str(size)


att_sweep(hp, scr_name, out_filename)
