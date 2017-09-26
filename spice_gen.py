class SpiceGenerator(object):
    def __init__(self, mesh_size):
        self.r_counter = 0
        self.v_counter = 0
        self.mesh_size = mesh_size

    def concat_name(self, name):
        if name != '':
            return '_'+name
        else:
            return ''

    def add_r(self, grid_idx1, grid_idx2, r, name=''):
        print("R"+str(self.r_counter)+self.concat_name(name),
              self.flatten_idx(grid_idx1),
              self.flatten_idx(grid_idx2),
              r)

    def add_v(self, grid_idx1, grid_idx2, v, name=''):
        if v > 0:
            print("V"+str(self.v_counter)+self.concat_name(name),
                  self.flatten_idx(grid_idx1),
                  self.flatten_idx(grid_idx2),
                  "PWL(0, "+str(v)+")")
        elif v < 0:
            print("V"+str(self.v_counter)+self.concat_name(name),
                  self.flatten_idx(grid_idx2),
                  self.flatten_idx(grid_idx1),
                  "PWL(0, "+str(-v)+")")

    def flatten_idx(self, idx):
        return idx[0]*self.mesh_size+idx[1]

    def add_block_comment(self, comment):
        self.add_newline()
        self.add_comment("")
        self.add_comment(comment)
        self.add_comment("")
        self.add_newline()

    def add_newline(self):
        print()

    def add_comment(self, comment):
        print("* "+comment)
