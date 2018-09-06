import numpy as np
import itertools as it
import sys
import math
import csv
import time

def write_difference_csv(mesh_size,N):
    with open('data/Analytical_Samples'+ str(mesh_size)+ 'Problem'+ str(N) + '.csv', 'r') as analytical, open('data/Spice_Voltage_Mesh'+ str(mesh_size)+ 'Problem'+ str(N) + '.csv', 'r') as spice, open('data/Difference_Mesh'+str(mesh_size)+'Problem'+str(N)+'.csv','w') as difference:
        analyticalReader = csv.reader(analytical)
        spiceReader = csv.reader(spice)
        differenceWriter = csv.writer(difference, delimiter = ',', quotechar='|', quoting = csv.QUOTE_MINIMAL)
        for rowA, rowS in zip(analyticalReader,spiceReader):
            tempRow=[]
            for A,S in zip(rowA,rowS):
                tempRow.append(abs(float(A) - float(S)))
            differenceWriter.writerow(tempRow)

def write_average_difference_csv(initial_mesh_size,mesh_size,N,start_time):
    with open('data/Difference_Mesh'+str(mesh_size)+'Problem'+str(N)+'.csv','r') as difference, open('data/AvgDifference_Mesh'+str(initial_mesh_size)+'ThroughMesh'+str(N)+'Problem'+str(N)+'.csv','a') as avgDifference:
        differenceReader = csv.reader(difference)
        avgDifferenceWriter = csv.writer(avgDifference, delimiter = ',', quotechar='|', quoting = csv.QUOTE_MINIMAL)
        rowCount = 0
        rowsTotalAverageDifference = 0
        meshAverageDifference=[]
        for row in differenceReader:
            rowTotalValue = 0
            numRowValues = 0
            for value in row:
                rowTotalValue += float(value)
                numRowValues += 1
            rowsTotalAverageDifference += rowTotalValue/numRowValues
            rowCount +=1
        gridAverageDifference = rowsTotalAverageDifference/rowCount
        meshAverageDifference.append(mesh_size)
        meshAverageDifference.append(gridAverageDifference)
        meshAverageDifference.append(time.time()-start_time)
        avgDifferenceWriter.writerow(meshAverageDifference)
