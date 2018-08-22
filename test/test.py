#! /usr/bin/env python
import sys
# print (sys.path)
import numpy as np
import itertools as it
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
# I am using python 3.6.1

# assume square grid
# N = 100
# source = (10, 10, 5, 5)
# sink = (55, 10, 5, 50)

#Not allowed, problem dimesions can only go to a minium of 10
#assume square problem dimensions of 5 by 5
#N = 5
#assume square problem dimensions of 25 by 25
#N = 25
N = 64
#assume square problem dimensions of 125 by 125
#N = 125

#Grid Layout
#y

#0
#1
#2
#3
#.
#.
#.
#24
#   0   1   2   ...  24 x


# Is the tuple laid onto the problem, or the mesh(grid) => The Problem, Good!
# Setup for the 25 problem, does not allow for variable source or sink => Bad
# Where is the source value set?
# I belive the source and sinks can not overlap => WRONG, Sources and Sinks CAN overlap
# I belive the source and sinks can not go out of bounds of the problem=> yes
# Index starts at 0 => Yes

# I am not going to overlap becuase I dont know how it effects values? =>
# The tuple is in form (left_pos(x axis), top_pos(y axis), width(x axis), height(y axis))
#top
#source = (0, 0, 25, 1)
source = (0, 0, 64, 1)
#left, bottom, right
#sink = [(0, 1, 1, 24),(0, 24, 24, 1),(24, 1, 1, 24)]
sink = [(0, 1, 1, 63),(0, 63, 63, 1),(63, 1, 1, 63)]
# sink = [(11, 2, 1, 10), (5, 15, 2, 2)]
# source = (4, 4, 4, 12)

cond_exp = -3
conductance = 10**cond_exp

#mesh dimension entered by user
mesh_size = int(sys.argv[1])
img_name = 'hm_{gr_sz}_{ms_sz}'

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
    write_function_csv(mesh_size,N)
    write_difference_csv(mesh_size,N)

    exponent+=1
    mesh_size += int(math.pow(2,exponent))

    if mesh_size > N:
        break

# eflow_data = energy_flow(m)
# print('Source in : {}'.format(eflow_data['src_in']))
# print('Source out: {}'.format(eflow_data['src_out']))
# print('Sink in   : {}'.format(eflow_data['snk_in']))
# print('Sink out  : {}'.format(eflow_data['snk_out']))

#print_current_table(m)
#write_current_table(m)
# print(N)
# print(mesh_size)
# name =  'Mesh'+ str(mesh_size)+ 'Problem:'+ str(N) + '.csv'
# print(name)
# with open('Mesh'+ str(mesh_size)+ 'Problem'+ str(N) + '.csv', 'w', newline = '') as csvfile:
#     writer = csv.writer(csvfile, delimiter = ',', quotechar='|', quoting = csv.QUOTE_MINIMAL)
#     for row in range(1,mesh_size):
#         row_currents = [l.ammeter.current for l  in m.links if is_in_row(l,row)]
#         writer.writerow(row_currents)

# write_current_csv(mesh_size,N,m)
# write_function_csv(mesh_size,N)
# write_difference_csv(mesh_size,N)

# with open('Analytical_Samples15Problem25.csv', 'r') as analytical, open('Spice_Mesh15Problem25.csv', 'r') as spice, open('Difference_Mesh15Problem25.csv','w') as difference:
#     analyticalReader = csv.reader(analytical)
#     spiceReader = csv.reader(spice)
#     differenceWriter = csv.writer(difference, delimiter = ',', quotechar='|', quoting = csv.QUOTE_MINIMAL)
#     for rowA, rowS in zip(analyticalReader,spiceReader):
#         tempRow=[]
#         for A,S in zip(rowA,rowS):
#             #print(abs(float(A) - float(S)))
#             tempRow.append(abs(float(A) - float(S)))
#         differenceWriter.writerow(tempRow)
        #print("good job")
    #print(difference)
        #file1 = t1.readlines()
        #file2 = t2.readlines()

# analytical = pandas.read_csv('Analytical_Samples15Problem25.csv')
# spice = pandas.read_csv('Spice_Mesh15Problem25.csv')
# difference = analytical - spice
# print

# with open('update.csv', 'w') as outFile:
#     for line in file1, file2:
#         outFile.write(abs(file1.line-file2.line))

# a=25
# b=25
# def f0(x,y):
#     return( x + y )
# def f(x,y):
#     return( (f0(x,y)/math.sinh((math.pi * b)/a)) * math.sin((math.pi * x / a))
#         + math.sinh(math.pi * y / a))

# row = 5
# fifth_row_currents = [l.ammeter.current for l in m.links if is_in_row(l,row)]
# print(fifth_row_currents)
# with open('data.csv', 'w', newline='') as csvfile:
#     spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
#     spamwriter.writerow(['Spam'] * 5 + ['Baked Beans'])
#     spamwriter.writerow(['Spam', 'Lovely Spam', 'Wonderful Spam'])
#     spamwriter.writerow(fifth_row_currents)


# current flows currently doesn't work correctly with the new grid
# logic. It would probably require some post processing on the grid
# itself (which sucks)

# if not surface3d:
#     plot_heatmap(m)
# else:
#     plot_surface(m)
