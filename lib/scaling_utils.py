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
from difference_utils import (write_difference_csv)

def write_exponential_scaling_csv(initial_mesh_size,N,source,sink,cond_exp,conductance):
    mesh_size = initial_mesh_size
    #powTwo = 2
    exponent = 0
    #RANGE (starting value, max value, iterator)
    while True:
        print("mesh_size: " + str(mesh_size))
        print("exponent: " + str(exponent))
    #for mesh_size in range(int(sys.argv[1]),N,int(math.pow(2,exponent))):
        #print(exponent)
        #exponent += 1
        #powerTwo = math.pow(2,exponent+=1)
        #print(powerTwo)
        hp = HeatProblem(N, source, sink, conductance, src_val=10.0)
        surface3d = False  # if True 3d surface plot is generated
        m = ROCModel(mesh_size)
        m.load_problem(hp)
        m.run_spice_solver()
        write_current_csv(mesh_size,N,m)
        # Missing Function transformationOFCurrentValues()
        write_function_csv(mesh_size,N)
        write_difference_csv(mesh_size,N)
        with open('data/Difference_Mesh'+str(mesh_size)+'Problem'+str(N)+'.csv','r') as difference, open('data/AvgDifference_Mesh'+str(initial_mesh_size)+'ThroughMesh'+str(N)+'Problem'+str(N)+'.csv','a') as avgDifference:
            differenceReader = csv.reader(difference)
            avgDifferenceWriter = csv.writer(avgDifference, delimiter = ',', quotechar='|', quoting = csv.QUOTE_MINIMAL)
            #for rowA, rowS in zip(analyticalReader,spiceReader):
            rowCount = 0
            rowsTotalAverageDifference = 0
            meshAverageDifference=[]
            for row in differenceReader:
                rowTotalValue = 0
                numRowValues = 0
                for value in row:
                    #print(abs(float(A) - float(S)))
                    #tempRow.append(abs(float(A) - float(S)))
                    rowTotalValue += float(value)
                    numRowValues += 1
                rowsTotalAverageDifference += rowTotalValue/numRowValues
                rowCount +=1
            gridAverageDifference = rowsTotalAverageDifference/rowCount
            meshAverageDifference.append(mesh_size)
            meshAverageDifference.append(gridAverageDifference)
            avgDifferenceWriter.writerow(meshAverageDifference)
        exponent+=1
        mesh_size += int(math.pow(2,exponent))
        if mesh_size > N:
            break
