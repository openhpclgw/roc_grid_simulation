import numpy as np
import itertools as it
import sys
import math
import csv

def write_function_csv(mesh_size,N):
    a=b=N
    with open('Analytical_Samples'+ str(mesh_size)+ 'Problem'+ str(N) + '.csv', 'w', newline = '') as csvfile:
        writer = csv.writer(csvfile, delimiter = ',', quotechar='|', quoting = csv.QUOTE_MINIMAL)
        #x = row & y = column
        for y in np.linspace(0,N,mesh_size,endpoint=False):
            currentRow=[]
            for x in np.linspace(0,N,mesh_size,endpoint=False):
                currentRow.append(f(x,y,a,b))
            writer.writerow(currentRow)

def f0(x,y):
    return( x + y )
def f(x,y,a,b):
    return( (f0(x,y)/math.sinh((math.pi * b)/a)) * math.sin((math.pi * x / a))
        + math.sinh(math.pi * y / a))
