import numpy as np
import itertools as it
import sys
import math
import csv

def write_difference_csv(mesh_size,N):
    with open('data/Analytical_Samples'+ str(mesh_size)+ 'Problem'+ str(N) + '.csv', 'r') as analytical, open('data/Spice_Mesh'+ str(mesh_size)+ 'Problem'+ str(N) + '.csv', 'r') as spice, open('data/Difference_Mesh'+str(mesh_size)+'Problem'+str(N)+'.csv','w') as difference:
        analyticalReader = csv.reader(analytical)
        spiceReader = csv.reader(spice)
        differenceWriter = csv.writer(difference, delimiter = ',', quotechar='|', quoting = csv.QUOTE_MINIMAL)
        for rowA, rowS in zip(analyticalReader,spiceReader):
            tempRow=[]
            for A,S in zip(rowA,rowS):
                #print(abs(float(A) - float(S)))
                tempRow.append(abs(float(A) - float(S)))
            differenceWriter.writerow(tempRow)
