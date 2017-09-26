class SpiceGenerator(object):

    def __init__(self, mesh_size):
        self.r_counter = 0
        self.v_counter = 0
        self.mesh_size = mesh_size
        self.__commentfrmt = '* {c}'
        self.__bcommentfrmt = '\n*\n* {c}\n*\n'
        self.__rfrmt = 'R{i}{uname} {n1} {n2} {r}'
        self.__vfrmt = 'V{i}{uname} {n1} {n2} PWL(0, {v})'

    #
    # codegen Functions
    #
    def add_r(self, grid_idx1, grid_idx2, r, name=''):
        print(self.__rfrmt.format(i=self.r_counter,
                                  uname=self.concat_name(name),
                                  n1=self.flatten_idx(grid_idx1),
                                  n2=self.flatten_idx(grid_idx2),
                                  r=r))
        self.r_counter += 1

    def add_v(self, grid_idx1, grid_idx2, v, name=''):
        if v > 0:
            print(self.__vfrmt.format(i=self.v_counter,
                                      uname=self.concat_name(name),
                                      n1=self.flatten_idx(grid_idx1),
                                      n2=self.flatten_idx(grid_idx2),
                                      v=v))
        elif v < 0:
            print(self.__vfrmt.format(i=self.v_counter,
                                      uname=self.concat_name(name),
                                      n1=self.flatten_idx(grid_idx2),
                                      n2=self.flatten_idx(grid_idx1),
                                      v=-v))
        self.v_counter += 1

    def add_block_comment(self, comment):
        print(self.__bcommentfrmt.format(c=comment))

    def add_comment(self, comment):
        print(self.__commentfrmt.format(c=comment))

    #
    # Utility Functions
    #
    def flatten_idx(self, idx):
        return idx[0]*self.mesh_size+idx[1]

    def concat_name(self, name):
        if name != '':
            return '_'+name
        else:
            return ''
