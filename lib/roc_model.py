import numpy as np
import itertools as it
from spice_gen import SpiceGenerator
from common import *


# Resistance is a simple spice object
class Resistance(object):
    def __init__(self, r, node1, node2, cntrs):
        self.uid = cntrs.r
        self.r = r
        self.node1 = node1
        self.node2 = node2
        # self.uname = uname
        cntrs.r += 1

# CurrentSource is a simple spice object
class CurrentSource(object):
    def __init__(self, i, node1, node2, cntrs):
        self.uid = cntrs.i
        self.i = i
        self.node1 = node1
        self.node2 = node2
        # self.uname = uname
        cntrs.i += 1

# BoundaryCond objects can be connected to any kind of node
class BoundaryCond(object):
    def __init__(self, v, node1, node2, cntrs): 
        self.uid = cntrs.v
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
        elif isinstance(node2, int): # Ground?
            if node2 == 0:
                self.node2 = '0'
        else:
            print('Unrecognized node' + str(node2))

        # self.uname = uname

        self.sname = ''

        cntrs.v += 1

class CurrentMeter(object):
    def __init__(self, node1, node2, cntrs):
        self.uid = cntrs.v # this counter logic needs to go
        self.node1 = node1
        self.node2 = node2
        self.current = 0.

        # Ammeters are used to inject potential within the mesh to do
        # current split analysis.
        self.v = 0

        cntrs.v += 1

    def set_bias(self, v):
        self.v = v

    def reset_bias(self):
        self.v = 0


# MeshResistances must connect two node block objects
class MeshResistance(object):
    def __init__(self, r, nodeblock1, nodeblock2, cntrs, norton=False): 
        # check if nodeblock2 comes "after" nodeblock1
        if nodeblock1 > nodeblock2:
            self.nodeblock1 = nodeblock1
            self.nodeblock2 = nodeblock2

        elif nodeblock2 > nodeblock1:
            self.nodeblock1 = nodeblock2
            self.nodeblock2 = nodeblock1
        else:
            print('Legs of MeshResistance cannot be same')

        # get corresponding terminals between NodeBlocks
        ts = nodeblock1.get_facing_subnodes(nodeblock2)
        if ts is None:
            print('MeshResistance can be created only between\
                    neighboring NodeBlocks')


        # if not norton:
        self.node1, self.node2 = ts
        # else:
            # set offset nodes
            # offset_node_frmt='{node}_off'
            # self.node1 = offset_node_frmt.format(node=ts[0])
            # self.node2 = offset_node_frmt.format(node=ts[1])

        self.orientation = 'H' if nodeblock1.is_h(nodeblock2) else 'V'
        self.resistance = Resistance(r, self.node1, self.node2,
                                     cntrs)

        # add internal subnode to attach the curmeter
        self.innodeid = str(self.resistance.uid)+'sub'

        # add curmeter
        self.curmeter = CurrentMeter(self.node1, self.innodeid, cntrs)
        self.resistance.node1 = self.innodeid

    def get_current(self):
        return self.curmeter.current

    def components(self):
        yield self.resistance
        yield self.curmeter

    def cur_direction(self, val):
        if val == 0:
            return 0
        if val < 0:
            return 'W' if self.orientation == 'H' else 'N'
        if val > 0:
            return 'E' if self.orientation == 'H' else 'S'


class NodeBlock(object):
    def __init__(self, coord, mesh_size, cntrs):
        def gen_node_name(coord, d=''):
            return 'N{}_{}{}'.format(coord[0], coord[1], d)

        self.coord = coord
        self.nodename = gen_node_name(coord)
        self.sname = self.nodename
        self.subnodes = {}
        self.curmeters = {}

        self.mesh_size = mesh_size

        self.potential = 0.

        # initial condition to be used in virtualized runs
        self.ic = 0.

        # TODO this should go to some sort of a global config file
        self.directions = ('E', 'W', 'N', 'S')

        for d in self.inward_directions():
            tmp_subnode = gen_node_name(coord, d)
            self.subnodes[d] = tmp_subnode
            self.curmeters[d] = CurrentMeter(node1=tmp_subnode,
                                       node2=self.nodename,
                                       cntrs=cntrs)

    def components(self):
        for d, a in self.curmeters.items():
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
        for d, a in self.curmeters.items():
            # print('\t', self.coord, ' ', d, ' ', a.current)
            if a.current > 0:
                ret += a.current
        return abs(ret)

    def sum_reduce_out_curs(self):
        ret = 0.
        for _, a in self.curmeters.items():
            if a.current < 0:
                ret += a.current
        return abs(ret)

    def aggregate_current_vector(self):
        a = self.curmeters
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

    def __repr__(self):
        return 'NodeBlock'+str(self.coord)

# TODO can we move the counter logic to the generators?
class CounterSet(object):
    def __init__(self):
        self.r = 0
        self.v = 0
        self.i = 0

class NortonLoop(object):
    def __init__(self, coord, nw, ne, sw, se, cntrs,
                 sides=None, boundary_cond=None):

        self.coord = coord
        self.nw = nw
        self.ne = ne
        self.sw = sw
        self.se = se
        self.nodeblocks = (nw, ne, sw, se)

        self.nodeblock_outdirs =  ((     'W', 'N'     ),
                                   ('E',      'N'     ),
                                   (     'W',      'S'),
                                   ('E',           'S'))

        self.edges = {'N': (self.nw, self.ne),
                      'E': (self.ne, self.se),
                      'S': (self.se, self.sw),
                      'W': (self.sw, self.nw)}

        self._components = []
        self.gnd_nbs = set()
        self.short_nbs = set()

        self.finalize(sides, boundary_cond, cntrs)

        self.current = 0.
        self.fixed_current = False
        self.current_branch = None

    def get_loop_current(self):
        if self.fixed_current:
            # we will return self.current at the end of the method
            # anyways, nothing to do here
            pass
        elif self.current_branch is not None:
            self.current = self.current_branch.get_current()
        else:
            # for nb in self.nodeblocks:
                # for d, curmeter in nb.curmeters.items():
                    # c = curmeter.current
                    # print(self.coord, ' ', nb.coord, ' ', d, ' ', c)
            cur = 0.
            for nb, dirs in zip(self.nodeblocks, self.nodeblock_outdirs):
                for d in dirs:
                    if d in nb.curmeters:
                        c = nb.curmeters[d].current
                        print(self.coord, ' ', nb.coord, ' ', d, ' ', c)
                        if c > 0:
                            cur += c
            self.current = cur
        return self.current

    def sum_reduce_in_curs(self):
        ret = 0.
        for nb in self.nodeblocks:
            cur = nb.sum_reduce_in_curs()
            ret += cur
        return abs(ret)

    def sum_reduce_out_curs(self):
        ret = 0.
        for nb in self.nodeblocks:
            cur = nb.sum_reduce_out_curs()
            ret += cur
        return abs(ret)

    def finalize(self, sides, boundary_cond, cntrs):
        # add the components depending on the location of the loop
        # and the location of sources and sinks

        # options
        # 1. Loop is on the side and source: Add current source on the
        # outer edge of the loop, 0 resistances on the other edges 2.
        # 2. Loop is on the side and sink: Connect the dangling
        # nodeblock nodes on the sideto ground, add 0 resistances on
        # other edges
        # 3. Loop is not on the side: Add 0 resistances on all edges

        # options 1 and 2 are handled in finalize_edge
        if sides is not None:
            if boundary_cond is None:
                bound_conn_type = 'short'
            else:
                bound_conn_type = boundary_cond
                if boundary_cond == 'source':
                    self.current = 1.0  # FIXME
                    self.fixed_current = True
                elif boundary_cond == 'sink':
                    self.current = 0.0
                    self.fixed_current = True
                elif boundary_cond == 'short':
                    # finalize edge will set self.current_branch once
                    # the MeshResistance object is initialized. Nothing
                    # to do here
                    pass

            for d,e in self.edges.items():
                if d in sides:
                    self.finalize_edge(e, bound_conn_type, cntrs)
                    # add only one soruce per loop
                    if boundary_cond == 'source':
                        bound_conn_type='short'
                else:
                    self.finalize_edge(e, 'no_change', cntrs)

        elif side is None:
            assert boundary_cond is None
            for _,e in self.edges.items():
                self.finalize_edge(e, 'no_change', cntrs)

    def finalize_edge(self, edge, conn_type, cntrs):
        print(edge)
        nb1 = edge[0]
        nb2 = edge[1]
        facing_nodes = nb1.get_facing_subnodes(nb2);
        n1 = facing_nodes[0]
        n2 = facing_nodes[1]
        oft = '{node}_off'
        if conn_type == 'no_change':
            pass

            # self.components.append(Resistance(r=0,
                                              # node1=n1,
                                              # node2=oft.format(node=n1),
                                              # cntrs=cntrs))
            # self.components.append(Resistance(r=0,
                                              # node1=n2,
                                              # node2=oft.format(node=n2),
                                              # cntrs=cntrs))
        elif conn_type == 'short':
            mesh_r = MeshResistance(r=0,
                                    nodeblock1=nb1,
                                    nodeblock2=nb2,
                                    cntrs=cntrs)
            self._components.append(mesh_r)
            self.short_nbs.add(nb1)
            self.short_nbs.add(nb2)
            self.current_branch = mesh_r
        elif conn_type == 'source':
            # self.components.append(BoundaryCond(v=0.001,
                                              # node1=n1,
                                              # node2=n2,
                                              # cntrs=cntrs))
            self._components.append(CurrentSource(i=1,
                                               node1=n1,
                                               node2=n2,
                                               cntrs=cntrs))
        elif conn_type == 'sink':
            self.gnd_nbs.add(nb1)
            self.gnd_nbs.add(nb2)
            # self.components.append(BoundaryCond(v=0,
                                              # node1=n1,
                                              # node2=0,
                                              # cntrs=cntrs))
            # self.components.append(BoundaryCond(v=0,
                                              # node1=n2,
                                              # node2=0,
                                              # cntrs=cntrs))
        else:
            print('Unrecognized conn_type: ' + conn_type)
            assert False

    def components(self):
        for c in self._components:
            if isinstance(c, MeshResistance):
                for cc in c.components():
                    yield cc
            else:
                yield c
            

class ROCModel(object):
    def __init__(self, N, norton=False):
        self.norton = norton
        self.h = N
        self.w = N

        if self.norton:
            self.h += 1
            self.w += 1

        self.r_counter = 0
        self.v_counter = 0
        self.mesh = np.zeros((self.meaningful_size(),
                              self.meaningful_size()))
        print(len(self.mesh))

        self.cntrs = CounterSet()

        self.nodes = []
        self.links = []
        if norton:
            self.loops = np.ndarray(shape=(self.w-1, self.w-1),
                                    dtype=NortonLoop)

        self.src_bboxs = []
        self.src_idxs = set()
        self.src = []

        self.snk_bboxs = []
        self.snk_idxs = set()
        self.snk = []

    def meaningful_size(self):
        if self.norton:
            return self.w-1
        else:
            return self.w

    def create_mesh(self, grid):

        def direct_copy(grid):
            self.mesh = grid.copy()
            print(len(self.mesh))
            self.exp_factor = 1.

        def interpolate_copy(grid):
            grid_size = len(grid)  # number of cells
            mesh_size = self.meaningful_size()  # number of nodes/loops

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
        print(hr)
        print(wr)

        self.nodes = [[NodeBlock(Coord(i, j), self.h, self.cntrs)
                      for j in wr] for i in hr]

    def init_links(self):
        # assert h = w

        # In norton circuit, we leave the edges un connected.
        # init_norton_loops creates all the loops, whose initializer
        # closes thos open circuits appropriately. This logic probably
        # needs to change if we start considering source/sink in the
        # middle meshes
        # if self.norton:
            # full_range = range(1, self.w-1True[
            # short_range = range(1, self.w-2)
        # else:
        full_range = range(self.w)
        short_range = range(self.w-1)

        # generate row resistors
        for i, j in it.product(full_range, short_range):
            if self.norton:
                if i == 0 or i == self.w-1:
                    continue
            self.links.append(MeshResistance(self.prob_conductance,
                                             self.nodes[i][j],
                                             self.nodes[i][j+1],
                                             self.cntrs,
                                             self.norton))

        # generate column resistors
        for i, j in it.product(short_range, full_range):
            if self.norton:
                if j == 0 or j == self.w-1:
                    continue
            self.links.append(MeshResistance(self.prob_conductance,
                                             self.nodes[i][j],
                                             self.nodes[i+1][j],
                                             self.cntrs,
                                             self.norton))

    def load_problem(self, hp):
        self.hp = hp
        grid = hp.gen_matrix()
        self.prob_conductance = hp.conductance
        self.init_mesh(grid)

    def get_snk_neighbor(my_idx):
        i = my_idx[0]
        j = my_idx[1]

        snk_nbors = set()

        if (i+0,j+1) in self.snk_idxs: snk_nbors.add('E')
        if (i+0,j-1) in self.snk_idxs: snk_nbors.add('W')
        if (i+1,j+0) in self.snk_idxs: snk_nbors.add('S')
        if (i-1,j+0) in self.snk_idxs: snk_nbors.add('N')

        return snk_nbors

    def iter_loops(self):
        assert self.norton

        full_range = range(self.w-1)  # mesh size is always node-wise
        for i, j in it.product(full_range, full_range):
            yield self.loops[i,j]


    # this needs to be called after source and sink is initialized
    def init_norton_loops(self):
        full_range = range(self.w-1)  # mesh size is always node-wise

        # SPICE doesn't like the same node to be grounded multiple
        # times. Once NortonLoop objects are created we query for the
        # nodes they need grounded. These nodes are stored in this set.
        # Elements must be NodeBlock objects 
        sink_nodes = set()
        shorted_nodes = set()

        for i, j in it.product(full_range, full_range):
            sides = set()
            if i == 0: sides.add('N')
            if i == self.w-2 : sides.add('S')
            if j == 0: sides.add('W')
            if j == self.w-2 : sides.add('E')

            if (i,j) in self.src_idxs:
                boundary_cond = 'source'
                print((i,j))
            elif (i,j) in self.snk_idxs:
                boundary_cond = 'sink'
            else:
                boundary_cond = None


            loop = NortonLoop(Coord(i,j),self.nodes[i+0][j+0],
                                         self.nodes[i+0][j+1],
                                         self.nodes[i+1][j+0],
                                         self.nodes[i+1][j+1],
                                         self.cntrs,
                                         sides,
                                         boundary_cond)

            shorted_nodes |= loop.short_nbs


            self.loops[i,j] = loop

        # for l in np.nditer(self.loops, flags=['refs_ok']):
        for i, j in it.product(full_range, full_range):
            l = self.loops[i,j]
            # now we have sink_nodes set storing all the nodeblocks that
            # needs to be grounded. But the question is which nodeblocks
            # that this loop owns are new (i.e. not been grounded
            # before)
            for nb in l.gnd_nbs - sink_nodes - shorted_nodes:
                l._components.append(BoundaryCond(v=0,
                                                    node1=nb,
                                                    node2=0,
                                                    cntrs=self.cntrs))

            # or the new nodeblocks into the nodes that have already
            # been grounded
            sink_nodes |= l.gnd_nbs

    # def calc_norton_loop_currents(self):
        # table = np.zeros((size, size))

        # iter_size = self.w-1
        # for l in self.loops:
            # if l

        

    def init_mesh(self, grid):
        self.create_mesh(grid)
        self.init_nodes()
        self.init_links()
        self.init_source(self.hp) # FIXME hp can go away
        self.init_sink(self.hp)
        if self.norton:
            self.init_norton_loops()
        else:
            self.init_ics(self.hp)  # TODO this can go away

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
            for s in hp.sources:
                tmp = mesh_bbox&s
                if tmp is not None:
                    self.src_bboxs.append(
                            tmp.lo(-vslice.left).to(-vslice.top))
            # print(self.src_bboxs)

        for s in self.src_bboxs:
            self.src_idxs |= {idx for idx in s}

        for i, j in self.src_idxs:
            self.src.append(BoundaryCond(self.mesh[i][j],
                                          self.nodes[i][j],
                                          0,
                                          self.cntrs))

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
                        self.src.append(BoundaryCond(self.mesh[i][j],
                                                      self.nodes[i][j],
                                                      0,
                                                      self.cntrs))

    # sources must be initialized before this is called
    # OR at least you need to make sure there is no overlap between
    # sources and sinks
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
            for s in hp.sinks:
                tmp = mesh_bbox&s
                if tmp is not None:
                    self.snk_bboxs.append(tmp.lo(-vslice.left).to(-vslice.top))

        for s in self.snk_bboxs:
            self.snk_idxs |= {idx for idx in s}

        self.snk_idxs -= self.src_idxs

        for i, j in self.snk_idxs:
            self.snk.append(BoundaryCond(v=0,
                                          node1=self.nodes[i][j],
                                          node2=0,
                                          cntrs=self.cntrs))

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

    def init_from_cache(self, filename):
        sg = SpiceGenerator(cache_only=True, filename=filename)
        sg.create_script(self)
        sg.get_results(self, cached_file=True)
        self.final_grid = self.create_grid()

    def run_spice_solver(self, filename='', 
                         cleanup=False, virtualize=False,
                         vstep_size=0):
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

        sg = SpiceGenerator(filename=filename)
        sg.create_script(self)
        sg.run()
        sg.get_results(self)
        grid = self.create_grid()
        count = 0
        if mesh_size < grid_size and virtualize:
            for s in grid_slices(vstep_size=vstep_size):
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

