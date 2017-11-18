import re
import roc_model as rm
import itertools as it
import numpy as np
import subprocess as sp


class SpiceGenerator(object):

    def __init__(self, filename=''):
        # init component counters
        self.r_counter = 0
        self.v_counter = 0
        self.tmp_folder = 'tmp'

        # create file handle
        if filename == '':
            filename = '__autogen_tmp'
        self.rel_in_path = '{}/{}.cir'.format(self.tmp_folder, filename)
        self.rel_out_path = '{}/{}.out'.format(self.tmp_folder, filename)

    def rm_tmp_files(self):
        sp.run('rm -f {} {}'.format(self.rel_in_path,
                                    self.rel_out_path),
               shell=True)

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
        self.file = open(self.rel_in_path, 'w')
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
        sp.run('ngspice -b {} -o {} >/dev/null'.format(self.rel_in_path,
                                            self.rel_out_path),
                                            shell=True)

    def get_results(self, roc_model):
        result_dict = self.generate_result_dict()

        for mr in roc_model.links:
            mr.ammeter.current = result_dict[mr.ammeter.sname]

        for node in roc_model.iter_nodes():
            for _,a in node.ammeters.items():
                a.current = result_dict[a.sname]

            sym = 'V('+node.sname+')'
            node.potential = result_dict[sym]

    def generate_result_dict(self):
        header_regexp = '^Index\s+time\s+(.*)$'
        value_regexp = '^0\s+\d[.]\d+e[+-]\d+\s+(\-?\d[.]\d+e[+-]\d+)'
        header_regobj = re.compile(header_regexp)
        value_regobj = re.compile(value_regexp)

        result_dict = {}

        looking_for_header = True
        name = ''
        val = 0.
        with open(self.rel_out_path) as f:
            for line in f:
                if looking_for_header:
                    m = header_regobj.match(line)
                    if m:
                        name = m.group(1).split('#')[0].strip().upper()
                        looking_for_header = False;
                else:
                    m = value_regobj.match(line)
                    if m:
                        val = float(m.group(1))
                        result_dict[name] = val
                        looking_for_header = True;

        return result_dict

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


