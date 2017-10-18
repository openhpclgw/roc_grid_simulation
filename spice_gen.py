from sys import stdout
import itertools as it
import numpy as np
import subprocess


class SpiceGenerator(object):

    def __init__(self, filename=''):
        # init component counters
        self.r_counter = 0
        self.v_counter = 0

        # create file handle
        if filename == '':
            self.file = stdout
        else:
            self.file = open('{}.cir'.format(filename), 'w')

    def set_id_pads(self, width):
        ps = str(width)

        # parametrized component id and node id rules
        self.cid_gr = '{i:0>'+ps+'}'  # set component id rule
        self.nid_gr = lambda i:\
            'N{n'+str(i)+'[0]:}_{n'+str(i)+'[1]:}{d'+str(i)+'}'

        # convenience names
        cg = self.cid_gr
        ng = self.nid_gr

        # define necessary spice grammar
        self.__commentfrmt = '* {c}'
        self.__bcommentfrmt = '\n*\n* {c}\n*'

        self.__nfrmt = 'N{n[0]:}_{n[1]:}'
        self.__rfrmt = 'R'+cg+'{uname} '+ng(1)+' '+ng(2)+' {r}'
        self.__vfrmt = 'V'+cg+'{uname} '+ng(1)+' '+ng(2)+' DC {v}'
        self.__pvfrmt = 'V'+cg+'{uname} N{n[0]:}_{n[1]:} 0 DC {v}'
        self.__tranfrmt = '.TRAN 1NS 11NS 10NS 10NS'
        self.__printfrmt = '.PRINT TRAN {typ}({symbol})'

    def create_script(self, mesh, conductance, hs):
        # mesh size
        self.mesh_size = len(mesh)
        self.set_id_pads(int(np.log10(self.mesh_size**2)+1))

        full_range = range(self.mesh_size)
        short_range = range(self.mesh_size-1)

        # generate row resistors
        self.add_block_comment("Row Resistors")
        for i, j in it.product(full_range, short_range):
            self.add_r((i, j), (i, j+1), conductance, horizontal=True)

        # generate column resistors
        self.add_block_comment("Column Resistors")
        for i, j in it.product(short_range, full_range):
            self.add_r((i, j), (i+1, j), conductance, horizontal=False)

        # I wanted to avoid adding external voltages that goes outside
        # the boundary and doesn't attach to anything. But spice seems
        # to go haywire if I do that
        add_extra_ammeters = True
        self.ammeters = [[[] for i in full_range] for j in full_range]
        # generate nodes with 4 0V sources
        self.add_block_comment("Node subcircuits")
        for i, j in it.product(full_range, full_range):
            self.add_comment("Node " + str((i, j)))
            if add_extra_ammeters or j > 0:
                self.ammeters[i][j].append(
                        self.add_v((i, j), (i, j), v=0, dir1='E'))
            if add_extra_ammeters or j < self.mesh_size-1:
                self.ammeters[i][j].append(
                        self.add_v((i, j), (i, j), v=0, dir1='W'))
            if add_extra_ammeters or i > 0:
                self.ammeters[i][j].append(
                        self.add_v((i, j), (i, j), v=0, dir1='N'))
            if add_extra_ammeters or i < self.mesh_size-1:
                self.ammeters[i][j].append(
                        self.add_v((i, j), (i, j), v=0, dir1='S'))

        # generate row voltages
        self.add_block_comment("Voltage Sources")
        for i, j in it.product(full_range, full_range):
            v = mesh[i][j]
            if v > 0:
                self.add_point_v((i, j), v)

        # generate heat sink
        for i,j in it.product(range(hs[1], hs[1]+hs[3]),
                              range(hs[0], hs[0]+hs[2])):
            self.add_point_v((i,j), mesh[i][j])

        # self.generate measurement/analysis components
        self.add_block_comment("Analysis code")
        self.add_transtmt()
        for i, j in it.product(full_range, full_range):
            for a in self.ammeters[i][j]:
                self.add_iprintstmt(a)
            self.add_vprintstmt((i,j))

        self.file.close()

    #
    # codegen Functions
    #
    def add_transtmt(self):
        self.gen(self.__tranfrmt.format())

    def add_iprintstmt(self, symbol):
        self.gen(self.__printfrmt.format(typ='I',
                                         symbol=symbol))

    def add_vprintstmt(self, n):
        self.gen(self.__printfrmt.format(typ='V',
                                         symbol=self.__nfrmt.format(n=n)))

    def add_r(self, grid_idx1, grid_idx2, r, horizontal, name=''):
        if horizontal:
            d1, d2 = 'E', 'W'
        else:
            d1, d2 = 'N', 'S'
        self.gen(self.__rfrmt.format(i=self.r_counter,
                                     uname=self.concat_name(name),
                                     n1=grid_idx1,
                                     d1=d1,
                                     n2=grid_idx2,
                                     d2=d2,
                                     r=r))

        self.r_counter += 1

    def add_v(self, grid_idx1, grid_idx2, v, dir1='', dir2='', name=''):
        ret = ''
        if v >= 0:
            ret = self.gen(
                    self.__vfrmt.format(i=self.v_counter,
                                        uname=self.concat_name(name),
                                        n1=grid_idx1,
                                        d1=dir1,
                                        n2=grid_idx2,
                                        d2=dir2,
                                        v=v))
        elif v < 0:
            ret = self.gen(
                    self.__vfrmt.format(i=self.v_counter,
                                        uname=self.concat_name(name),
                                        n1=grid_idx2,
                                        d1=dir1,
                                        n2=grid_idx1,
                                        d2=dir2,
                                        v=-v))
        self.v_counter += 1
        return ret

    def add_point_v(self, grid_idx, v, name=''):
        self.gen(self.__pvfrmt.format(i=self.v_counter,
                                      uname=self.concat_name(name),
                                      n=grid_idx,
                                      v=v))
        self.v_counter += 1

    def add_block_comment(self, comment):
        self.gen(self.__bcommentfrmt.format(c=comment))

    def add_comment(self, comment):
        self.gen(self.__commentfrmt.format(c=comment))

    #
    # simulator utils
    #
    def run(self):
        # FIXME set file names/hierarchy better
        subprocess.run('ngspice -b test.cir -o test.out', shell=True)

    def get_results(self, hsrc, hsnk):
        w = self.mesh_size
        h = self.mesh_size
        results = np.zeros((h, w))
        # FIXME parse in a more pythonic way
        grep_cmd = 'grep -i {sym} test.out -A2 \
                    | tail -n1 \
                    | tr "\t" " " \
                    | cut -d" " -f"3"'
        vid = 0
        U = np.zeros((self.mesh_size, self.mesh_size))
        V = np.zeros((self.mesh_size, self.mesh_size))
        ein = np.zeros((self.mesh_size, self.mesh_size))
        eout = np.zeros((self.mesh_size, self.mesh_size))
        for i, j in it.product(range(h), range(w)):
            def accumulate_energy(val):
                if val < 0:
                    eout[i][j] += val
                else:
                    ein[i][j] += val
            a = self.ammeters[i][j]
            east = float(subprocess.check_output(
                                   grep_cmd.format(sym=a[0]),
                                   shell=True))
            accumulate_energy(east)
            U[i][j] += -east
            
            west = float(subprocess.check_output(
                                   grep_cmd.format(sym=a[1]),
                                   shell=True))
            accumulate_energy(west)
            U[i][j] += west

            north = float(subprocess.check_output(
                                   grep_cmd.format(sym=a[2]),
                                   shell=True))
            accumulate_energy(north)
            V[i][j] += -north
            
            south = float(subprocess.check_output(
                                   grep_cmd.format(sym=a[3]),
                                   shell=True))
            accumulate_energy(south)
            V[i][j] += south

            sym = 'v\('+self.__nfrmt.format(n=(i,j))+'\)'
            tmp_val = float(subprocess.check_output(
                                       grep_cmd.format(sym=sym),
                                       shell=True))

            if tmp_val > 0:
                results[i][j] += tmp_val
            vid += 1

        sum_out = 0.
        for i,j in it.product(range(hsrc[1], hsrc[1]+hsrc[3]),
                              range(hsrc[0], hsrc[0]+hsrc[2])):
            sum_out += eout[i][j]

        sum_in = 0.
        for i,j in it.product(range(hsnk[1], hsnk[1]+hsnk[3]),
                              range(hsnk[0], hsnk[0]+hsnk[2])):
            sum_in+= ein[i][j]

        return results, U, V, sum_out, sum_in

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
        return s.split()[0]
