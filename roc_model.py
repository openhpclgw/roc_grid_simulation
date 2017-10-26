import numpy as np
import itertools as it
from spice_gen import SpiceGenerator

# class NodeName(object):
    # def __init__(self, name):
        # self.name=name

class Coord(object):
    def __init__(self, i, j):
        self.i = i
        self.j = j
        self.__coord = (i,j)

    def __getitem__(self, key):
        return self.__coord[key]

    def distance(self, other):
        return abs(self.i-other.i) + abs(self.j-other.j)

    def is_neighbor(self, other):
        return self.distance(other) == 1

# Resistance is a simple spice object
class Resistance(object):
    counter = 0
    def __init__(self, r, node1, node2, uname=''):
        self.uid = Resistance.counter
        self.r = r
        self.node1 = node1
        self.node2 = node2
        self.uname = uname
        Resistance.counter += 1

# VoltageSource objects can be connected to any kind of node
class VoltageSource(object):
    counter = 0
    def __init__(self, v, node1, node2, uname=''):
        self.uid = VoltageSource.counter
        self.v = v
        self.node1 = node1
        self.node2 = node2
        self.uname = uname

        self.sname = ''

        VoltageSource.counter += 1



class Ground(VoltageSource):
    def __init__(self, nid):
        VoltageSource.__init__(self, v=0, node1=node, node2=0)

class Ammeter(VoltageSource):
    def __init__(self, node1, node2):
        VoltageSource.__init__(self, v=0, node1=node1, node2=node2)

# MeshResistances must connect two node block objects
class MeshResistance(object):
    def __init__(self, r, nodeblock1, nodeblock2, uname=''):
        # self.r = r
        self.nodeblock1 = nodeblock1
        self.nodeblock2 = nodeblock2
        # self.uname = uname

        ts = nodeblock1.get_facing_subnodes(nodeblock2)
        self.node1, self.node2 = ts

        self.resistance = Resistance(r, self.node1, self.node2)

        # add internal subnode to attach the ammeter
        self.innodeid = str(self.resistance.uid)+'sub'

        # add ammeter
        self.ammeter = Ammeter(self.node1, self.innodeid)
        self.resistance.node1 = self.innodeid

    def components(self):
        yield self.resistance
        yield self.ammeter

class NodeBlock(object):
    def __init__(self, coord):
        def gen_node_name(coord, d=''):
            return 'N{}_{}{}'.format(coord[0], coord[1], d)

        self.coord = coord
        self.nodename = gen_node_name(coord)
        self.subnodes = {}
        self.ammeters = {}

        # TODO this should go to some sort of a global config file
        self.directions = ('E', 'W', 'N', 'S')

        for d in self.directions:
            tmp_subnode = gen_node_name(coord, d)
            self.subnodes[d] = tmp_subnode
            self.ammeters[d] = Ammeter(tmp_subnode, self.nodename)

    def components(self):
        for d,a in self.ammeters.items():
            yield a

    def get_facing_subnodes(self, other):
        if not self.is_neighbor(other):
            return None

        if self.coord.i == other.coord.i:
            if self.coord.j < other.coord.j:
                return (self.subnodes['E'], other.subnodes['W'])
            else:
                return (self.subnodes['W'], other.subnodes['E'])
        else:
            if self.coord.i < other.coord.i:
                return (self.subnodes['S'], other.subnodes['N'])
            else:
                return (self.subnodes['N'], other.subnodes['S'])
        
    def is_neighbor(self, other):
        return self.coord.is_neighbor(other.coord)



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
            return 1.

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
            return extrp_factor
        #
        # end of interpolate
        #

        # grid_h = len(grid)
        grid_w = len(grid[0])

        if grid_w <= self.w:
            return direct_copy(grid)
        else:  # problem is larger than the mesh
            return interpolate_copy(grid)
    #
    # end of create_mesh
    #

    def init_nodes(self):
        hr = range(self.h)
        wr = range(self.w)
        self.nodes = [[NodeBlock(Coord(i,j)) for j in wr] for i in hr]


    def init_links(self, conductance):
        # assert h = w
        full_range = range(self.w)
        short_range = range(self.w-1)

        # generate row resistors
        self.links = []
        for i, j in it.product(full_range, short_range):
            self.links.append(MeshResistance(conductance,
                                             self.nodes[i][j],
                                             self.nodes[i][j+1]))

        # generate column resistors
        for i, j in it.product(short_range, full_range):
            self.links.append(MeshResistance(conductance,
                                             self.nodes[i][j],
                                             self.nodes[i+1][j]))


    def load_problem(self, grid, conductance):
        self.exp_factor = self.create_mesh(grid)
        self.init_nodes()
        self.init_links(conductance)
        self.prob_conductance = conductance

    def run_spice_solver(self, hp):
        self.load_problem(hp.gen_matrix(), hp.conductance)
        # create the heatsink zone
        # print(self.exp_factor)
        extrp_hsnk = (int(hp.sink[0]/self.exp_factor),
                      int(hp.sink[1]/self.exp_factor),
                      int(np.ceil(hp.sink[2]/self.exp_factor)),
                      int(np.ceil(hp.sink[3]/self.exp_factor)))
        extrp_hsrc = (int(hp.source[0]/self.exp_factor),
                      int(hp.source[1]/self.exp_factor),
                      int(np.ceil(hp.source[2]/self.exp_factor)),
                      int(np.ceil(hp.source[3]/self.exp_factor)))
        print(extrp_hsrc)
        print(extrp_hsnk)
        sg = SpiceGenerator('test')
        sg.create_script(self.mesh, self.prob_conductance, extrp_hsnk,
                self)  # FIXME
        sg.run()
        return sg.get_results(extrp_hsrc, extrp_hsnk, self)
