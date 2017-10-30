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

    def __str__(self):
        return str(self.__coord)

    def distance(self, other):
        return abs(self.i-other.i) + abs(self.j-other.j)

    def is_neighbor(self, other):
        return self.distance(other) == 1

    def is_h(self, other):
        return abs(self.j-other.j)==1

    # this should be enough to run the logic for now
    def __gt__(self, other):
        return self.i>other.i or (self.i==other.i and self.j>other.j)

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
    def __init__(self, v, node1, node2='0', uname=''):
        self.uid = VoltageSource.counter
        self.v = v
        if isinstance(node1, str):
            self.node1 = node1
        elif isinstance(node1, NodeBlock):
            self.node1 = node1.nodename
        else:
            print('Unrecognized node' + str(node1))

        if isinstance(node2, str):
            self.node2 = node2
        elif isinstance(node2, NodeBlock):
            self.node2 = node2.nodename
        else:
            print('Unrecognized node' + str(node2))

        self.uname = uname

        self.sname = ''

        VoltageSource.counter += 1



class Ground(VoltageSource):
    def __init__(self, node):
        VoltageSource.__init__(self, v=0, node1=node, node2='0')

class Ammeter(VoltageSource):
    def __init__(self, node1, node2):
        VoltageSource.__init__(self, v=0, node1=node1, node2=node2)
        self.current = 0.

    def set_bias(self, v):
        self.v = v

    def reset_bias(self):
        self.v = 0

# MeshResistances must connect two node block objects
class MeshResistance(object):
    def __init__(self, r, nodeblock1, nodeblock2, uname=''):
        # check if nodeblock2 comes "after" nodeblock1
        if not (nodeblock2> nodeblock1):
            print('2nd NodeBlock of MeshResistance must come after the\
                  1st in ij ordering')

        # get corresponding terminals between NodeBlocks
        ts = nodeblock1.get_facing_subnodes(nodeblock2)
        if ts == None:
            print('MeshResistance can be created only between\
                    neighboring NodeBlocks')

        # self.r = r
        self.nodeblock1 = nodeblock1
        self.nodeblock2 = nodeblock2

        self.node1, self.node2 = ts

        self.orientation = 'H' if nodeblock1.is_h(nodeblock2) else 'V'
        self.resistance = Resistance(r, self.node1, self.node2, uname)

        # add internal subnode to attach the ammeter
        self.innodeid = str(self.resistance.uid)+'sub'

        # add ammeter
        self.ammeter = Ammeter(self.node1, self.innodeid)
        self.resistance.node1 = self.innodeid

    def components(self):
        yield self.resistance
        yield self.ammeter
    
    def cur_direction(self, val):
        if val == 0:
            return 0
        if val < 0:
            return 'W' if self.orientation=='H' else 'N'
        if val > 0:
            return 'E' if self.orientation=='H' else 'S'


class NodeBlock(object):
    def __init__(self, coord):
        def gen_node_name(coord, d=''):
            return 'N{}_{}{}'.format(coord[0], coord[1], d)

        self.coord = coord
        self.nodename = gen_node_name(coord)
        self.sname = self.nodename
        self.subnodes = {}
        self.ammeters = {}

        self.potential = 0.

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

    def is_h(self, other):
        return self.coord.is_h(other.coord)

    def __gt__(self, other):
        return self.coord>other.coord

    def sum_reduce_in_curs(self):
        ret = 0.
        for _,a in self.ammeters.items():
            if a.current > 0:
                ret += a.current
        return abs(ret)

    def sum_reduce_out_curs(self):
        ret = 0.
        for _,a in self.ammeters.items():
            if a.current < 0:
                ret += a.current
        return abs(ret)

    def aggregate_current_vector(self):
        u = self.ammeters['W'].current-self.ammeters['E'].current
        v = self.ammeters['S'].current-self.ammeters['N'].current
        return (u,v)

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
            # print(grid)
            # print(self.mesh)
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

    # for convenience. I think this can be refactored out if nodes was a
    # numpy array
    def iter_nodes(self):
        for i,j in it.product(range(self.h), range(self.w)):
            yield self.nodes[i][j]

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


    def load_problem(self, hp):
        grid = hp.gen_matrix()
        conductance = hp.conductance
        self.exp_factor = self.create_mesh(grid)
        self.init_nodes()
        self.init_links(conductance)
        self.init_source(hp)
        self.init_sink(hp)
        self.prob_conductance = conductance

    def init_source(self, hp):
        # what to do in case of indivisible sizes?
        self.src_bboxs = [(int(s[0]/self.exp_factor),
                         int(s[1]/self.exp_factor),
                         int(np.ceil(s[2]/self.exp_factor)),
                         int(np.ceil(s[3]/self.exp_factor)))
                         for s in hp.source_iter()]
        self.src_idxs = set()
        self.src = []
        for s in self.src_bboxs:
            self.src_idxs |= {idx for idx in self.__iter_bbox(s)}

        for i,j in self.src_idxs:
            self.src.append(VoltageSource(self.mesh[i][j],
                                          self.nodes[i][j]))

    def init_sink(self, hp):
        # what to do in case of indivisible sizes?
        self.snk_bboxs = [(int(s[0]/self.exp_factor),
                         int(s[1]/self.exp_factor),
                         int(np.ceil(s[2]/self.exp_factor)),
                         int(np.ceil(s[3]/self.exp_factor))) 
                         for s in hp.sink_iter()]
        self.snk_idxs = set()
        self.snk = []
        for s in self.snk_bboxs:
            self.snk_idxs |= {idx for idx in self.__iter_bbox(s)}

        for i,j in self.snk_idxs:
            self.snk.append(Ground(self.nodes[i][j]))

    def __iter_bbox(self, bbox):
        for i,j in it.product(range(bbox[1], bbox[1]+bbox[3]),
                              range(bbox[0], bbox[0]+bbox[2])):
            yield (i,j)

    def src_nodeblocks(self):
        for i,j in self.src_idxs:
            yield self.nodes[i][j]

    def snk_nodeblocks(self):
        for i,j in self.snk_idxs:
            yield self.nodes[i][j]

    def run_spice_solver(self, hp, cleanup=False):
        sg = SpiceGenerator()
        sg.create_script(self)
        sg.run()
        sg.get_results(self)
        if cleanup:
            sg.rm_tmp_files()



