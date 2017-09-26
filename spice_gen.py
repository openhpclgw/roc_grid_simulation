from sys import stdout
import itertools as it
import numpy as np


class SpiceGenerator(object):

    def __init__(self, filename=''):
        # define necessary spice grammar
        self.__commentfrmt = '* {c}'
        self.__bcommentfrmt = '\n*\n* {c}\n*'

        # these two are set once we know how much to pad
        self.__rfrmt = ''
        self.__vfrmt = ''

        # init component counters
        self.r_counter = 0
        self.v_counter = 0

        # create file handle
        if filename == '':
            self.file = stdout
        else:
            self.file = open('{}.cir'.format(filename), 'w')

    def __call__(self, mesh, conductance):
        # mesh size
        self.mesh_size = len(mesh)
        self.set_id_pads(int(np.log10(self.mesh_size**2)+1))

        full_range = range(self.mesh_size)
        short_range = range(self.mesh_size-1)

        # generate row resistors
        self.add_block_comment("Row Resistors")
        for i, j in it.product(full_range, short_range):
            self.add_r((i, j), (i, j+1), conductance)

        # generate column resistors
        self.add_block_comment("Column Resistors")
        for i, j in it.product(short_range, full_range):
            self.add_r((i, j), (i+1, j), conductance)

        # generate row voltages
        self.add_block_comment("Row Voltage Sources")
        for i, j in it.product(full_range, short_range):
            self.add_v((i, j), (i, j+1), (mesh[i][j]-mesh[i][j+1]))

        # generate row voltages
        self.add_block_comment("Column Voltage Sources")
        for i, j in it.product(short_range, full_range):
            self.add_v((i, j), (i+1, j), (mesh[i][j]-mesh[i+1][j]))

        # self.generate measurement/analysis components

    #
    # codegen Functions
    #
    def add_r(self, grid_idx1, grid_idx2, r, name=''):
        self.gen(self.__rfrmt.format(i=self.r_counter,
                                     uname=self.concat_name(name),
                                     n1=self.flatten_idx(grid_idx1),
                                     n2=self.flatten_idx(grid_idx2),
                                     r=r))
        self.r_counter += 1

    def add_v(self, grid_idx1, grid_idx2, v, name=''):
        if v > 0:
            self.gen(self.__vfrmt.format(i=self.v_counter,
                                         uname=self.concat_name(name),
                                         n1=self.flatten_idx(grid_idx1),
                                         n2=self.flatten_idx(grid_idx2),
                                         v=v))
        elif v < 0:
            self.gen(self.__vfrmt.format(i=self.v_counter,
                                         uname=self.concat_name(name),
                                         n1=self.flatten_idx(grid_idx2),
                                         n2=self.flatten_idx(grid_idx1),
                                         v=-v))
        self.v_counter += 1

    def add_block_comment(self, comment):
        self.gen(self.__bcommentfrmt.format(c=comment))

    def add_comment(self, comment):
        self.gen(self.__commentfrmt.format(c=comment))

    #
    # Utility Functions
    #
    def flatten_idx(self, idx):
        return idx[0]*self.mesh_size+idx[1]

    def concat_name(self, name):
        if name != '':
            return '_{}'.format(name)
        else:
            return ''

    def gen(self, s):
        self.file.write('{}\n'.format(s))

    def set_id_pads(self, width):
        ps = str(width)
        self.__rfrmt = 'R{i:0>'+ps+'}{uname} {n1:>'+ps+'} {n2:>'+ps+'} {r}'
        self.__vfrmt = 'V{i:0>'+ps+'}{uname} {n1:>'+ps+'} {n2:>'+ps+'} PWL(0, {v})'
