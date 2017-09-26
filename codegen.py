class SpiceModel(object):
    def __init__(self, grid, cond_table):
        self.grid = grid
        self.cond_table = cond_table
        self.h= len(self.grid)
        self.w = len(self.grid[0])
        self.r_counter = 0
        self.v_counter = 0

    def to_spice(self):

        # generate row resistors
        self.add_block_comment("Row Resistors")
        for i in range(self.w):
            self.add_comment("Row " + str(i) + " resistors")
            for j in range(self.w-1):
                self.add_r((i,j), (i,j+1), 1, str(self.r_counter))
                self.r_counter += 1

        # generate column resistors
        self.add_block_comment("Column Resistors")
        for i in range(self.w-1):
            self.add_comment("Column " + str(i) + " resistors")
            for j in range(self.w):
                self.add_r((i,j), (i+1,j), 1, str(self.r_counter))
                self.r_counter += 1

        # generate row voltages
        self.add_block_comment("Row Voltage Sources")
        for i in range(self.w):
            self.add_comment("Row " + str(i) + " voltage sources")
            for j in range(self.w-1):
                self.add_v((i,j), (i,j+1),
                           (self.grid[i][j]-self.grid[i][j+1]),
                           str(self.v_counter))
                self.v_counter += 1

        # generate row voltages
        self.add_block_comment("Column Voltage Sources")
        for i in range(self.w-1):
            self.add_comment("Row " + str(i) + " voltage sources")
            for j in range(self.w):
                self.add_v((i,j), (i+1,j),
                           (self.grid[i][j]-self.grid[i+1][j]),
                           str(self.v_counter))
                self.v_counter += 1

        # generate measurement/analysis components


    def add_r(self, grid_idx1, grid_idx2, r, name):
        print("R"+name,
              self.flatten_idx(grid_idx1),
              self.flatten_idx(grid_idx2),
              r)

    def add_v(self, grid_idx1, grid_idx2, v, name):
        if v > 0:
            print("V"+name,
                    self.flatten_idx(grid_idx1),
                    self.flatten_idx(grid_idx2),
                    "PWL(0, "+str(v)+")")
        elif v < 0:
            print("V"+name,
                    self.flatten_idx(grid_idx2),
                    self.flatten_idx(grid_idx1),
                    "PWL(0, "+str(-v)+")")


    def flatten_idx(self, idx):
        return idx[0]*self.w+idx[1]

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

