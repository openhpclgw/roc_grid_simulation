from heat_problem import HeatProblem
from att_sweep_helper import *

def resistance(size, n1, n2):
    return 0.001*multiplier

size = mesh_size
source = (int(size/2)-1, int(size/2)-1, 2, 2)
sink = [(0, 0, 1, size),  # left
        (size-1, 0, 1, size),  # right
        (1, 0, size-2, 1),      # top
        (1, size-1, size-2, 1)]   # bottom
hp = HeatProblem(size, source, sink, resistance)

working_dir = 'test/attenuation_sweep/case2/size'+str(mesh_size)+'/'
scr_name = working_dir+'opt_case2_'+str(size)+'_att'+str(base_att)
out_filename = working_dir+'att'+str(base_att)+'/{}_'+str(size)

att_sweep(hp, scr_name, out_filename, working_dir)
