import csv
import math
import numpy as np


    
def parse_comsol_csv(filename):
    with open(filename, newline='\n') as csvfile:
        reading_data = False
        comsol_reader = csv.reader(csvfile, delimiter=',')
        offset = 0.
        for r in comsol_reader:
            if not reading_data:
                if r[0] == '% Nodes':
                    num_nodes = int(r[1])
                    tmp_grid_size = math.sqrt(num_nodes)
                    grid_size = int(tmp_grid_size)
                    if grid_size != tmp_grid_size:
                        print('Grid is not square')
                    grid = np.zeros((grid_size, grid_size))

                elif r[0] == '% x':
                    reading_data = True
            else: # now we are reading the data
                def to_ij(xy):
                    return (int(grid_size*(xy[1]-offset)),
                            int(grid_size*(xy[0]-offset)))

                tmp_x = float(r[0])
                tmp_y = float(r[1])
                tmp_val = float(r[2])
                if offset == 0.:
                    offset = tmp_x

                grid[to_ij((tmp_x, tmp_y))] = tmp_val

        return grid
                    
print(parse_comsol_csv('test.csv'))
