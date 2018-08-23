import numpy as np
import itertools as it
import sys
import math
import csv
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import (aggregate_current_vectors,
                            print_current_table,
                            write_current_csv,
                            is_in_row,
                            energy_flow,
                            plot_surface,
                            plot_heatmap)
import csv
import math
from analytical_solutions import (write_function_csv,f0,f)
import pandas
from difference_utils import (write_difference_csv,write_average_difference_csv)

def write_exponential_scaling_csv(initial_mesh_size,N,source,sink,cond_exp,conductance):
    mesh_size = initial_mesh_size
    exponent = 0
    while True:
        print("mesh_size: " + str(mesh_size))
        print("exponent: " + str(exponent))
        hp = HeatProblem(N, source, sink, conductance, src_val=10.0)
        surface3d = False  # if True 3d surface plot is generated
        m = ROCModel(mesh_size)
        m.load_problem(hp)
        m.run_spice_solver()
        write_current_csv(mesh_size,N,m)
        #!!! Missing Function transformationOFCurrentValues() !!!
        write_function_csv(mesh_size,N)
        write_difference_csv(mesh_size,N)
        write_average_difference_csv(initial_mesh_size,mesh_size,N)
        exponent+=1
        mesh_size += int(math.pow(2,exponent))
        if mesh_size > N:
            break

def write_linear_scaling_csv(initial_mesh_size,N,source,sink,cond_exp,conductance):
    mesh_size = initial_mesh_size
    while True:
        print("mesh_size: " + str(mesh_size))
        #print("exponent: " + str(exponent))
        hp = HeatProblem(N, source, sink, conductance, src_val=10.0)
        surface3d = False  # if True 3d surface plot is generated
        m = ROCModel(mesh_size)
        m.load_problem(hp)
        m.run_spice_solver()
        write_current_csv(mesh_size,N,m)
        #!!! Missing Function transformationOFCurrentValues() !!!
        write_function_csv(mesh_size,N)
        write_difference_csv(mesh_size,N)
        write_average_difference_csv(initial_mesh_size,mesh_size,N)
        mesh_size += 1
        if mesh_size > N:
            break
