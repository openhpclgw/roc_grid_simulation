import numpy as np
import itertools as it
from spice_gen import SpiceGenerator


class ROCModel(object):
    def __init__(self, N):
        self.h = N
        self.w = N
        self.r_counter = 0
        self.v_counter = 0
        self.mesh = np.zeros((self.h, self.w))

    def create_mesh(self, grid):

        def direct_copy(grid):
            self.mesh = grid.copy()

        def interpolate_copy(grid):
            grid_size = len(grid)  # number of cells
            mesh_size = self.w  # number of nodes

            # establish grid indices and coords
            grid_x = np.arange(0, grid_size, 1)
            grid_y = np.arange(0, grid_size, 1)
            # grid_x_coords = grid_x + .5
            # grid_y_coords = grid_x + .5

            # ditto for mesh
            # mesh_x = np.arange(0, mesh_size, 1)
            # mesh_y = np.arange(0, mesh_size, 1)

            # calculate necessary scaling/offset constants
            extrp_factor = float(grid_size)/float(mesh_size)
            phys_offset = extrp_factor/2.

            # establish cell bounds
            g_bnds = np.arange(0, grid_size+1, 1)
            m_bnds = np.arange(0, grid_size+1, extrp_factor)

            # for a given cell return it's bounds
            def cell_to_bound(i, j, bounds):
                return (bounds[j+1], bounds[i], bounds[j], bounds[i+1])

            # intesect two bounds objects. any negativity means
            # rectangles don't have an intersection
            def intersect(b1, b2):
                return (min(b1[0], b2[0]), max(b1[1], b2[1]),
                        max(b1[2], b2[2]), min(b1[3], b2[3]))

            # calcute area of an intersection
            def area(b):
                if b[0]-b[2] < 0 or b[3]-b[1] < 0:
                    return 0.
                return (b[0]-b[2])*(b[3]-b[1])

            # iterate nearest mesh cell indices for a given grid index
            def nearest_mesh_cells(g_i, g_j):

                def __nmc(g_io, g_jo):
                    yield (int(np.floor(g_io/extrp_factor))-1,
                           int(np.floor(g_jo/extrp_factor))-1)
                    yield (int(np.floor(g_io/extrp_factor))-1,
                           int(np.floor(g_jo/extrp_factor)))
                    yield (int(np.floor(g_io/extrp_factor)),
                           int(np.floor(g_jo/extrp_factor))-1)
                    yield (int(np.floor(g_io/extrp_factor)),
                           int(np.floor(g_jo/extrp_factor)))

                g_io = g_i+phys_offset
                g_jo = g_j+phys_offset

                for tmp_cell in __nmc(g_io, g_jo):
                    if tmp_cell[0] >= 0 and tmp_cell[0] < mesh_size:
                        if tmp_cell[1] >= 0 and tmp_cell[1] < mesh_size:
                            yield tmp_cell

            # iterate over the grid and calculate weighted sums
            # for g_i in grid_x:
                # for g_j in grid_y:
            for (g_i, g_j) in it.product(grid_x, grid_y):
                for (m_i, m_j) in nearest_mesh_cells(g_i, g_j):
                    self.mesh[m_i][m_j] += (grid[g_i][g_j]*area(intersect(
                                            cell_to_bound(g_i, g_j, g_bnds),
                                            cell_to_bound(m_i, m_j, m_bnds))))

            # average the mesh
            # self.mesh = [m/(extrp_factor**2) for m in self.mesh]
            self.mesh = self.mesh/(extrp_factor**2)
            # print(extrp_factor)
            print(grid)
            print(self.mesh)
        #
        # end of interpolate
        #

        # grid_h = len(grid)
        grid_w = len(grid[0])

        if grid_w <= self.w:
            direct_copy(grid)
        else:  # problem is larger than the mesh
            interpolate_copy(grid)
    #
    # end of create_mesh
    #

    def to_spice(self, grid, conductance):
        self.create_mesh(grid)
        sg = SpiceGenerator()
        sg(self.mesh, conductance)
