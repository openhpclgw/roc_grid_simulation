from sys import stdout
import roc_model as rm
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
        self.__v2frmt = 'V'+cg+'{uname} {nn1} {nn2} DC {v}'
        self.__r2frmt = 'R'+cg+'{uname} {nn1} {nn2} {r}'
        self.__pvfrmt = 'V'+cg+'{uname} N{n[0]:}_{n[1]:} 0 DC {v}'
        self.__tranfrmt = '.TRAN 1NS 501NS 100NS 100NS'
        self.__printfrmt = '.PRINT TRAN {typ}({symbol})'

    def create_script(self, mesh, conductance, hs, roc_model):
        # mesh size
        self.mesh_size = len(mesh)
        self.set_id_pads(int(np.log10(self.mesh_size**2*4)+1))

        full_range = range(self.mesh_size)
        short_range = range(self.mesh_size-1)

        # generate row resistors
        self.add_block_comment('Resistors')
        for r in roc_model.links:
            for c in r.components():
                self.component_codegen(c)

        # generate nodes 
        self.add_block_comment("Node subcircuits")
        for i, j in it.product(full_range, full_range):
            self.add_comment("Node " + str((i, j)))
            for c in roc_model.nodes[i][j].components():
                self.component_codegen(c)

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
        # for s in roc_model.snk_nodes():
            # self.component_codegen(c)

        # self.generate measurement/analysis components
        self.add_block_comment("Analysis code")
        self.add_transtmt()
        for i, j in it.product(full_range, full_range):
            for d,a in roc_model.nodes[i][j].ammeters.items():
                self.add_iprintstmt(a.sname)
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
            d1, d2 = 'S', 'N'
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

    def get_results(self, hsrc, hsnk, roc_model):
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
        rcurs = []
        for i, j in it.product(range(h), range(w)):
            def accumulate_energy(val):
                if val < 0:
                    eout[i][j] += val
                else:
                    ein[i][j] += val
            a = roc_model.nodes[i][j].ammeters
            # print([v.sname for k,v in a.items()])
            # print(a['E'].sname)
            # print(grep_cmd.format(sym=a['E'].sname))
            east = float(subprocess.check_output(
                                   grep_cmd.format(sym=a['E'].sname),
                                   shell=True))
            if j+1 < w:
                direction = 0 if east==0 else 'W' if east>0 else 'E'
                rcurs.append(((i,j),   # terminal 1
                              (i,j+1),  # terminal 2
                              abs(east),  # current
                              direction))  # +1: West, -1: East
            accumulate_energy(east)
            U[i][j] += -east
            
            west = float(subprocess.check_output(
                                   grep_cmd.format(sym=a['W'].sname),
                                   shell=True))
            accumulate_energy(west)
            U[i][j] += west

            north = float(subprocess.check_output(
                                   grep_cmd.format(sym=a['N'].sname),
                                   shell=True))
            accumulate_energy(north)
            V[i][j] += -north
            
            south = float(subprocess.check_output(
                                   grep_cmd.format(sym=a['S'].sname),
                                   shell=True))
            if i+1 < h:
                direction = 0 if south==0 else 'N' if south>0 else 'S'
                rcurs.append(((i,j),   # terminal 1
                              (i+1,j),  # terminal 2
                              abs(south),  # current
                              direction))
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
        print(hsrc)
        # for i,j in it.product(range(hsrc[1], hsrc[1]+hsrc[3]),
                              # range(hsrc[0], hsrc[0]+hsrc[2])):
            # print((i,j))
        for i,j in roc_model.src_nodes():
            sum_out += eout[i][j]

        sum_in = 0.
        for i,j in it.product(range(hsnk[1], hsnk[1]+hsnk[3]),
                              range(hsnk[0], hsnk[0]+hsnk[2])):
            print((i,j))
            sum_in+= ein[i][j]

        return results, U, V, sum_out, sum_in, rcurs

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

    def add_r2(self, r):
        ret = self.gen(self.__r2frmt.format(i=r.uid,
                                     uname=self.concat_name(r.uname),
                                     nn1=r.node1,
                                     nn2=r.node2,
                                     r=r.r))

        self.r_counter += 1

    def add_v2(self, v):
        ret = ''
        if v.v >= 0:
            ret = self.gen(
                    self.__v2frmt.format(i=v.uid,
                                        uname=self.concat_name(v.uname),
                                        nn1=v.node1,
                                        nn2=v.node2,
                                        v=v.v))
        elif v.v < 0:
            print("How?")
            self.__v2frmt.format(i=v.uid,
                                uname=self.concat_name(v.uname),
                                nn1=v.node2,
                                nn2=v.node1,
                                v=v.v)
        self.v_counter += 1
        return ret

    def component_codegen(self, c):
        if isinstance(c, rm.VoltageSource):
            tmp_name = self.add_v2(c)
        elif isinstance(c, rm.Resistance):
            tmp_name = self.add_r2(c)
        else:
            tmp_name = ''
            print("error")
        c.sname = tmp_name
        
    def gen(self, s):
        self.file.write('{}\n'.format(s))
        return s.split()[0]


