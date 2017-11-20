import numpy as np
import itertools as it
from spice_gen import SpiceGenerator
from common import *


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
        if not (nodeblock2 > nodeblock1):
            print('2nd NodeBlock of MeshResistance must come after the\
                  1st in ij ordering')

        # get corresponding terminals between NodeBlocks
        ts = nodeblock1.get_facing_subnodes(nodeblock2)
        if ts is None:
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
            return 'W' if self.orientation == 'H' else 'N'
        if val > 0:
            return 'E' if self.orientation == 'H' else 'S'


class NodeBlock(object):
    def __init__(self, coord, mesh_size):
        def gen_node_name(coord, d=''):
            return 'N{}_{}{}'.format(coord[0], coord[1], d)

        self.coord = coord
        self.nodename = gen_node_name(coord)
        self.sname = self.nodename
        self.subnodes = {}
        self.ammeters = {}

        self.mesh_size = mesh_size

        self.potential = 0.

        # initial condition to be used in virtualized runs
        self.ic = 0.

        # TODO this should go to some sort of a global config file
        self.directions = ('E', 'W', 'N', 'S')

        for d in self.inward_directions():
            tmp_subnode = gen_node_name(coord, d)
            self.subnodes[d] = tmp_subnode
            self.ammeters[d] = Ammeter(tmp_subnode, self.nodename)

    def components(self):
        for d, a in self.ammeters.items():
            yield a

    def inward_directions(self):
        if self.coord.i > 0:
            yield 'N'
        if self.coord.i < self.mesh_size-1:
            yield 'S'
        if self.coord.j > 0:
            yield 'W'
        if self.coord.j < self.mesh_size-1:
            yield 'E'

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
        return self.coord > other.coord

    def sum_reduce_in_curs(self):
        ret = 0.
        for _, a in self.ammeters.items():
            if a.current > 0:
                ret += a.current
        return abs(ret)

    def sum_reduce_out_curs(self):
        ret = 0.
        for _, a in self.ammeters.items():
            if a.current < 0:
                ret += a.current
        return abs(ret)

    def aggregate_current_vector(self):
        a = self.ammeters
        curs = {}
        for d in self.directions:
            if d in a:
                if a[d].current < 0:
                    curs[d] = abs(a[d].current)
                else:
                    curs[d] = 0.
            else:
                curs[d] = 0.
        u = curs['E']-curs['W']
        v = curs['N']-curs['S']
        return u, v


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
            self.exp_factor = 1.

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
            self.exp_factor = extrp_factor
            # print(extrp_factor)
            # print(grid)
            # print(self.mesh)
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
        for i, j in it.product(range(self.h), range(self.w)):
            yield self.nodes[i][j]

    def init_nodes(self):
        hr = range(self.h)
        wr = range(self.w)
        self.nodes = [[NodeBlock(Coord(i, j), self.h)
                      for j in wr] for i in hr]

    def init_links(self):
        # assert h = w
        full_range = range(self.w)
        short_range = range(self.w-1)

        # generate row resistors
        self.links = []
        for i, j in it.product(full_range, short_range):
            self.links.append(MeshResistance(self.prob_conductance,
                                             self.nodes[i][j],
                                             self.nodes[i][j+1]))

        # generate column resistors
        for i, j in it.product(short_range, full_range):
            self.links.append(MeshResistance(self.prob_conductance,
                                             self.nodes[i][j],
                                             self.nodes[i+1][j]))

    def load_problem(self, hp):
        self.hp = hp
        grid = hp.gen_matrix()
        self.prob_conductance = hp.conductance
        self.init_mesh(grid)

    def init_mesh(self, grid):
        self.create_mesh(grid)
        self.init_nodes()
        self.init_links()
        self.init_source(self.hp) # FIXME hp can go away
        self.init_sink(self.hp)
        self.init_ics(self.hp)

    def init_virtualized_mesh(self, grid, vslice):
        print('Initializing virtualized mesh')
        # TODO TODO TODO here the data from the previous grid is not
        # transferred to the new virtualized grid. Only sources and
        # sinks are -correctly- created but the node potentials from the
        # initial interpolated run is lost
        self.create_mesh(grid)
        self.init_nodes()
        self.init_links()
        self.init_source(self.hp, vslice) # FIXME hp can go away
        self.init_sink(self.hp, vslice)
        self.init_ics(self.hp, vslice)
        print('Initialized virtualized mesh')

    def create_grid(self):
        ef = self.exp_factor
        if int(ef) != ef:
            print('Warning: Extrapolation factor is not integer')

        interm_g_s = int(self.h*ef)
        interm_grid = np.zeros((interm_g_s, interm_g_s))
        for i,j in it.product(range(interm_g_s), range(interm_g_s)):
            interm_grid[i][j] = self.nodes[int(i/ef)][int(j/ef)].potential

        return interm_grid

    def clear_mesh(self):
        self.nodes.clear()
        self.links.clear()
        self.src.clear()
        self.src_bboxs.clear()
        self.src_idxs.clear()
        self.ics = np.zeros((self.w, self.w))
        self.snk.clear()
        self.snk_bboxs.clear()
        self.snk_idxs.clear()


    def init_source(self, hp, vslice=None):
        # what to do in case of indivisible sizes?
        if vslice == None:
            self.src_bboxs = [BoundingBox(int(s[0]/self.exp_factor),
                              int(s[1]/self.exp_factor),
                              int(np.ceil(s[2]/self.exp_factor)),
                              int(np.ceil(s[3]/self.exp_factor)))
                              for s in hp.sources]
        else:
            mesh_bbox = BoundingBox(vslice.left,vslice.top,self.w,self.h)
            self.src_bboxs = []
            for s in hp.sources:
                tmp = mesh_bbox&s
                if tmp is not None:
                    self.src_bboxs.append(
                            tmp.lo(-vslice.left).to(-vslice.top))
            print(self.src_bboxs)

        self.src_idxs = set()
        self.src = []
        for s in self.src_bboxs:
            self.src_idxs |= {idx for idx in s}

        for i, j in self.src_idxs:
            self.src.append(VoltageSource(self.mesh[i][j],
                                          self.nodes[i][j]))

    def init_ics(self, hp, vslice=None):
        ms = self.w
        for i,j in it.product(range(ms), range(ms)):
            if (i,j) not in self.src_idxs:
                self.nodes[i][j].ic = self.mesh[i][j]

        # if this was a vslice and there was no source within the slice
        # then use the boundaries as voltage source
        if vslice is not None:
            if len(self.src_idxs) == 0:
                for i, j in self.boundaries():
                    if (i,j) not in self.snk_idxs:
                        self.src.append(VoltageSource(self.mesh[i][j],
                                                      self.nodes[i][j]))

    def init_sink(self, hp, vslice=None):
        # what to do in case of indivisible sizes?
        if vslice==None:
            self.snk_bboxs = [BoundingBox(int(s[0]/self.exp_factor),
                              int(s[1]/self.exp_factor),
                              int(np.ceil(s[2]/self.exp_factor)),
                              int(np.ceil(s[3]/self.exp_factor)))
                              for s in hp.sinks]
        else:
            mesh_bbox = BoundingBox(vslice.left,vslice.top,self.w,self.h)
            self.snk_bboxs = []
            for s in hp.sinks:
                tmp = mesh_bbox&s
                if tmp is not None:
                    self.snk_bboxs.append(tmp.lo(-vslice.left).to(-vslice.top))

        self.snk_idxs = set()
        self.snk = []
        for s in self.snk_bboxs:
            self.snk_idxs |= {idx for idx in s}

        for i, j in self.snk_idxs:
            self.snk.append(Ground(self.nodes[i][j]))

    def src_nodeblocks(self):
        for i, j in self.src_idxs:
            yield self.nodes[i][j]

    def snk_nodeblocks(self):
        for i, j in self.snk_idxs:
            yield self.nodes[i][j]

    def node_potentials(self):
        ms = self.h
        return np.array([[self.nodes[j][i].potential
                        for i in range(ms)] for j in range(ms)])

    # iterate mesh boundaries wioth no particular order
    def boundaries(self):
        ms = self.h
        for i in range(ms):
            yield (i,0)
            yield (i,ms-1)

        for j in range(1, ms-1):
            yield (0,j)
            yield (ms-1,j)

    def run_spice_solver(self, cleanup=False, virtualize=False):
        mesh_size = self.h
        grid_size = self.hp.N

        def get_grid_slice(grid, bbox):
            grid_slice = np.zeros((bbox.width, bbox.height))

            i_off = bbox.top
            j_off = bbox.left

            for i,j in it.product(range(mesh_size), range(mesh_size)):
                grid_slice[i][j] = grid[i+i_off][j+j_off]

            return grid_slice

        def grid_slices(vstep_size=0):
            ep = self.exp_factor
            if vstep_size == 0:
                vstep_size = mesh_size
            slice_range = range(0, grid_size-mesh_size+1, vstep_size)
            count = 0
            for l,t in it.product(slice_range, slice_range):
                yield BoundingBox(l,t,mesh_size, mesh_size)
                count += 1

        def update_grid_slice(grid, data, bbox):
            i_off = bbox.top
            j_off = bbox.left

            for i,j in it.product(range(mesh_size), range(mesh_size)):
                grid[i+i_off][j+j_off] = data[i][j]

        def grid_debug():
            for i in range(grid_size):
                for j in range(grid_size):
                    print('{:>5.2f} '.format(grid[i][j]), end='')
                print()

        sg = SpiceGenerator()
        sg.create_script(self)
        sg.run()
        sg.get_results(self)
        grid = self.create_grid()
        count = 0
        if mesh_size < grid_size and virtualize:
            for s in grid_slices(vstep_size=1):
                # print('Slice')
                # print(s)
                # print('PRE')
                # grid_debug()
                self.clear_mesh()
                self.init_virtualized_mesh(get_grid_slice(grid, s), s)
                sg.create_script(self, count)
                sg.run(count)
                sg.get_results(self, count)
                # print('GENERATED')
                # print(self.node_potentials())
                update_grid_slice(grid, self.node_potentials(), s)
                # print('POST')
                # grid_debug()
                count += 1

        self.final_grid = grid

        if cleanup:
            sg.rm_tmp_files()

