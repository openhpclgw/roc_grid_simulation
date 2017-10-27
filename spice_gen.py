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

    def create_script(self, roc_model):
        # mesh size
        self.mesh_size = len(roc_model.mesh)
        self.set_id_pads(int(np.log10(self.mesh_size**2*4)+1))

        full_range = range(self.mesh_size)
        short_range = range(self.mesh_size-1)

        # generate resistors
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

        # generate source
        self.add_block_comment("Sources")
        for v in roc_model.src:
            self.component_codegen(v)

        # generate sink
        self.add_block_comment("Sinks")
        for s in roc_model.snk:
            self.component_codegen(s)


        # self.generate measurement/analysis components
        self.add_block_comment("Analysis code")
        self.add_transtmt()
        for i, j in it.product(full_range, full_range):
            for _,a in roc_model.nodes[i][j].ammeters.items():
                self.add_iprintstmt(a.sname)
            self.add_vprintstmt((i,j))

        for mr in roc_model.links:
            self.add_iprintstmt(mr.ammeter.sname)

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

    def get_results(self, roc_model):
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

        for mr in roc_model.links:
            cur = float(subprocess.check_output(
                                   grep_cmd.format(sym=mr.ammeter.sname),
                                   shell=True))
            rcurs.append((mr.nodeblock1.coord, mr.nodeblock2.coord,
                abs(cur), 
                mr.cur_direction(cur)))

        for i, j in it.product(range(h), range(w)):
            node = roc_model.nodes[i][j]
            a = node.ammeters
            east = float(subprocess.check_output(
                                   grep_cmd.format(sym=a['E'].sname),
                                   shell=True))
            node.accumulate_energy(east)
            U[i][j] += -east
            
            west = float(subprocess.check_output(
                                   grep_cmd.format(sym=a['W'].sname),
                                   shell=True))
            node.accumulate_energy(west)
            U[i][j] += west

            north = float(subprocess.check_output(
                                   grep_cmd.format(sym=a['N'].sname),
                                   shell=True))
            node.accumulate_energy(north)
            V[i][j] += -north
            
            south = float(subprocess.check_output(
                                   grep_cmd.format(sym=a['S'].sname),
                                   shell=True))
            node.accumulate_energy(south)
            V[i][j] += south

            sym = 'v\('+self.__nfrmt.format(n=(i,j))+'\)'
            tmp_val = float(subprocess.check_output(
                                       grep_cmd.format(sym=sym),
                                       shell=True))

            if tmp_val > 0:
                results[i][j] += tmp_val
            vid += 1

        sum_out = 0.
        for i,j in roc_model.src_nodes():
            sum_out += roc_model.nodes[i][j].eout

        sum_in = 0.
        for i,j in roc_model.snk_nodes():
            sum_in += roc_model.nodes[i][j].ein

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


