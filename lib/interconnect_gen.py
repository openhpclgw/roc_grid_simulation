import re
import roc_model as rm
import itertools as it
import numpy as np
import subprocess as sp

coord_to_sch_ratio = 10.
sch_offset = 0.8

class InterconnectGenerator(object):

    def __init__(self, cache_only=False, filename=''):
        # init component counters
        self.r_counter = 0
        self.v_counter = 0
        self.osc_counter = 0
        self.sparam_counter = 0
        self.tmp_folder = 'tmp'

        self.cache_only = cache_only

        # create file handle
        if filename == '':
            self.user_fname = False
            self.filename = '__autogen_tmp'
        else:
            self.user_fname = True
            self.filename = filename

        self.counters = { 'ring': 0, 'y': 0, 'pwm': 0, 'cwl': 0,
                          'fiber': 1000 }

    def rel_out_path(self, suffix=''):
        if self.user_fname:
            return '{}.out'.format(self.filename)
        else:
            return '{}/{}_{}.out'.format(self.tmp_folder,
                                         self.filename,
                                         suffix)

    def rel_in_path(self, suffix=''):
        if self.user_fname:
            return '{}.cir'.format(self.filename)
        else:
            return '{}/{}_{}.cir'.format(self.tmp_folder,
                                         self.filename,
                                         suffix)
    def rel_lsf_path(self, suffix=''):
        if self.user_fname:
            return '{}.lsf'.format(self.filename)
        else:
            return '{}/{}_{}.lsf'.format(self.tmp_folder,
                                         self.filename,
                                         suffix)

    def rm_tmp_files(self):
        sp.run('rm -f {} {}'.format(self.rel_in_path(),
                                    self.rel_out_path()),
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

        self.__schfrmt = ' sch_x={sch_x:4.2f} sch_y={sch_y:4.2f} sch_f=f lay_x=0 lay_y=0 {custom}'
        self.__nfrmt = 'N{n[0]:}_{n[1]:}'
        self.__rfrmt = 'R'+cg+'{uname} '+ng(1)+' '+ng(2)+' {r}'
        self.__vfrmt = 'V'+cg+'{uname} '+ng(1)+' '+ng(2)+' DC {v}'
        self.__v2frmt = 'X_CWL_'+cg+' {nn1} \"CW Laser\" \"internal seed\"=7434 '+self.__schfrmt
        self.__r2frmt = 'X_FIBER_'+cg+'{uname} {nn1} {nn2} \"Optical Linear Fiber\" attenuation={r} ' + self.__schfrmt
        self.__pvfrmt = 'V'+cg+'{uname} N{n[0]:}_{n[1]:} 0 DC {v}'

        self.__oscfrmt = 'X_OOSC_'+cg+' {nn2} \"Optical Oscilloscope\"  '+ self.__schfrmt
        self.__sparamfrmt = 'X_SPAR_'+cg+' '+ng(1)+' '+ng(2)+' '+ng(3)+' '+ng(4)+' \"Optical N Port S-Parameter\" \"s parameters filename\"=\"{coord}.txt\" '+ self.__schfrmt

        self.__tranfrmt = '.TRAN 1NS 501NS 100NS 100NS'
        self.__printfrmt = '.PRINT TRAN {typ}({symbol})'
        self.__icfrmt = '.IC V({node}) {val}'
        self.sparam_lsf_format = ('addelement("Optical N Port S-Parameter");\n'
                                  'set("name","SPAR_'+cg+'");\n'
                                  'set("load from file",true);\n'
                                  'set("s parameters filename","spar.txt");\n'
                                  'set("x position", {sch_x:4.2f});\n'
                                  'set("y position", {sch_y:4.2f});\n'
                                  '{conns}\n')
        self.conn_frmt = ('connect("SPAR_'+cg+'","{this_port}",'
                                  ' "{other}", "{other_port}");\n')
        self.__pwmfrmt = 'X_OPWM_'+cg+' {nn1} \"Optical Power Meter\"'+self.__schfrmt
        self.__yfrmt = 'X_Y_'+cg+' {nn1} {nn2} {nn3} \"Waveguide Y Branch\"'+self.__schfrmt
        self.__ringfrmt = 'X_RING_'+cg+' {nn1} {nn2} \"Single Bus Ring Resonator\" frequency=230T '+self.__schfrmt

    def create_script(self, roc_model, suffix=''):
        if not self.cache_only:
            self.file = open(self.rel_in_path(suffix), 'w')
            self.lsffile = open(self.rel_lsf_path(suffix), 'w')
        # mesh size
        self.mesh_size = len(roc_model.mesh)
        self.set_id_pads(int(np.log10(self.mesh_size**2*4)+1))

        self.add_block_comment('Auto-generated for Interconnect')

        self.gen_header()

        full_range = range(self.mesh_size)
        short_range = range(self.mesh_size-1)

        # generate resistors
        self.add_block_comment('Resistors')
        for r in roc_model.links:
            for c in r.components():
                self.component_codegen(c, r)

        # generate nodes 
        self.add_block_comment("Node subcircuits")
        for i, j in it.product(full_range, full_range):
            n = roc_model.nodes[i][j]
            self.add_comment("Node " + str((i, j)))
            for c in n.components():
                self.component_codegen(c, n, model=roc_model)

            # generate the initial condition
            if n.ic != 0.:
                self.ic_codegen(n.sname, n.ic)

        # generate source
        self.add_block_comment("Sources")
        for v in roc_model.src:
            self.component_codegen(v, model=roc_model)

        # generate sink
        self.add_block_comment("Sinks")
        for s in roc_model.snk:
            self.component_codegen(s, model=roc_model)


        # self.generate measurement/analysis components
        self.add_block_comment("Analysis code")
        # self.add_transtmt()
        # for i, j in it.product(full_range, full_range):
            # for _,a in roc_model.nodes[i][j].ammeters.items():
                # self.add_iprintstmt(a.sname)
            # self.add_vprintstmt((i,j))

        # for mr in roc_model.links:
            # self.add_iprintstmt(mr.ammeter.sname)

        if not self.cache_only:
            self.file.close()
            self.lsffile.close()

    def nodename_to_coord(self, name):
        coord_str = name[1:].split('_')
        return (int(coord_str[0]), int(coord_str[1]))

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
    def run(self, suffix=''):
        sp.run('ngspice -b {} -o {} >/dev/null'.format(
                                            self.rel_in_path(suffix),
                                            self.rel_out_path(suffix)),
                                            shell=True)

    def get_results(self, roc_model, suffix='',
                    cached_file=False):

        result_dict = self.generate_result_dict(suffix, cached_file)

        for mr in roc_model.links:
            mr.ammeter.current = result_dict[mr.ammeter.sname]

        for node in roc_model.iter_nodes():
            for _,a in node.ammeters.items():
                a.current = result_dict[a.sname]

            sym = 'V('+node.sname+')'
            node.potential = result_dict[sym]

    def generate_result_dict(self, suffix, cached_file=False):
        header_regexp = '^Index\s+time\s+(.*)$'
        value_regexp = '^0\s+\d[.]\d+e[+-]\d+\s+(\-?\d[.]\d+e[+-]\d+)'
        header_regobj = re.compile(header_regexp)
        value_regobj = re.compile(value_regexp)

        result_dict = {}

        looking_for_header = True
        name = ''
        val = 0.
        fname = self.filename+'.out' if cached_file else self.rel_out_path(suffix)

        with open(fname) as f:
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

    def midpoint(self, sch_coord1, sch_coord2):
        if sch_coord1[0] == sch_coord2[0]:
            return (sch_coord1[0], (sch_coord1[1]+sch_coord2[1])/2)
        else:
            return ((sch_coord1[0]+sch_coord2[0])/2, sch_coord1[1])

    def node_sch_coord(self, coord):
        return (coord[1]*coord_to_sch_ratio,
                coord[0]*coord_to_sch_ratio)

    def add_r2(self, r, parent=None):
        custom = ''
        if parent is not None:
            node1_sch = self.node_sch_coord(parent.nodeblock1.coord)
            node2_sch = self.node_sch_coord(parent.nodeblock2.coord)
            sch_x, sch_y = self.midpoint(node1_sch, node2_sch)

            if parent.orientation == 'V':
                custom = '"rotated"=true'
        else:
            sch_x, sch_y = 0, 0

        

        ret = self.gen(self.__r2frmt.format(i=r.uid,
                                     uname='', # FIXME
                                     nn1=r.node1,
                                     nn2=r.node2,
                                     r=r.r,
                                     sch_x=sch_x,
                                     sch_y=0-sch_y,
                                     custom=custom))

        self.r_counter += 1
        return ret


    # self.node_to_bc = {}
    def add_v2(self, v):
        ret = ''

        sch_x, sch_y = self.node_sch_coord(
                            self.nodename_to_coord(v.node1))
        
        if sch_x == 0:
            sch_x -= sch_offset*2
        else:
            sch_x += sch_offset*2
        if sch_y == 0:
            sch_y -= sch_offset*2
        else:
            sch_y += sch_offset*2



        if v.v >= 0:
            ret = self.gen(
                    self.__v2frmt.format(i=v.uid,
                                         uname='',
                                         nn1=v.node1,
                                         nn2=v.node2,
                                         v=v.v,
                                         sch_x=sch_x,
                                         sch_y=0-sch_y,
                                         custom=''))
        elif v.v < 0:
            print("How?")
            self.__v2frmt.format(i=v.uid,
                                uname='',
                                nn1=v.node2,
                                nn2=v.node1,
                                v=v.v)
        self.v_counter += 1

        
        # self.node_to_bc[v.node1] = ret
        return ret

    def add_osc(self, osc, parent=None):
        custom=''
        if parent is not None:
            if isinstance(parent, rm.MeshResistance):
                node1_sch = self.node_sch_coord(parent.nodeblock1.coord)
                node2_sch = self.node_sch_coord(parent.nodeblock2.coord)
                sch_x, sch_y = self.midpoint(node1_sch, node2_sch)

                if parent.orientation == 'H':
                    sch_y += sch_offset*0.7
                else:
                    custom = '"rotated"=true'
                    sch_x -= sch_offset*1.1

            elif isinstance(parent, rm.NodeBlock):
                print(parent.coord)
                sch_x, sch_y = self.node_sch_coord(parent.coord)
                d = osc.node1[-1]
                if d == 'E':
                    sch_x += sch_offset
                    sch_y += sch_offset
                elif d == 'W':
                    sch_x -= sch_offset
                    sch_y -= sch_offset
                elif d == 'N':
                    sch_x += sch_offset
                    sch_y -= sch_offset
                elif d == 'S':
                    sch_x -= sch_offset
                    sch_y += sch_offset

        else:
            sch_x, sch_y = 0, 0
        ret = ''
        ret = self.gen(
                self.__oscfrmt.format(i=osc.uid,
                                     uname='',
                                     nn1=osc.node1,
                                     nn2=osc.node2,
                                     sch_x=sch_x,
                                     sch_y=0-sch_y,
                                     custom=custom))
        self.osc_counter += 1
        return ret

    def add_sparam(self, conn, sch_x=0, sch_y=0):
        ret = ''
        ret = self.gen(
                self.__sparamfrmt.format(i=conn.uid,
                                     uname='',
                                     n1=conn.coord,
                                     n2=conn.coord,
                                     n3=conn.coord,
                                     n4=conn.coord,
                                     d1='E',
                                     d2='W',
                                     d3='N',
                                     d4='S',
                                     coord=str(conn.coord),
                                     sch_x=sch_x,
                                     sch_y=0-sch_y,
                                     custom=''))
        self.sparam_counter += 1
        return ret

    def gen_header(self):
        self.file.write('{}\n'.format( (
        '.MODEL "CW Laser" ellipticity=0 "enable RIN"=0\n'
        '+ "orthogonal identifier 1"=1 "label 1"="X" seed=1\n'
        '+ "automatic seed"=1 linewidth=0 "orthogonal identifier 2"=2\n'
        '+ power=0.001 frequency=193.1T "reference power"=0.001 RIN=-140\n'
        '+ phase=0 "label 2"="Y" azimuth=0\n'
        '+ "number of output signals"={%number of output signals%}\n'
        '+ "output signal mode"={%output signal mode%}\n'
        '+ "sample rate"={%sample rate%} "time window"={%time window%}\n'
        '\n'
        '.MODEL "Optical Linear Fiber" "filter fit tolerance"=0.001\n'
        '+ "number of fir taps"=64 "dispersion slope"=80\n'
        '+ "single tap filter"=0 "run diagnostic"=0\n'
        '+ configuration="bidirectional" "reference frequency"=193.1T\n'
        '+ "diagnostic size"=1024 length=0.001\n'
        '+ "estimate number of taps"=1 dispersion=16u modes="X,Y" \n'
        '+ "window function"="rectangular"\n'
        '\n'
        '.MODEL "Optical N Port S-Parameter" "passivity tolerance"=1u\n'
        '+ configuration="bidirectional" "digital filter type"="FIR"\n'
        '+ passivity="ignore" "number of iir taps"=4\n'
        '+ "filter fit tolerance"=0.05 "maximum number of iir taps"=20\n'
        '+ order=1 "number of fir taps"=64 "remove disconnected ports"=0\n'
        '+ "run diagnostic"=0 "filter fit rolloff"=0.05\n'
        '+ reciprocity="ignore" "window function"="rectangular"\n'
        '+ "filter fit number of iterations"=50\n'
        '+ "estimate number of taps"=0 "load from file"=1\n'
        '+ "diagnostic size"=1024 temperature={%temperature%}\n'
        '\n'
        '.MODEL "Optical Oscilloscope" "input signal selection"="last"\n'
        '+ "refresh length"=1024 seed=1 refresh=1 "automatic seed"=1\n'
        '+ "display memory length"=2048 "limit display memory"=1\n'
        '+ "power unit"="W" frequency=193.1T "plot format"="power"\n'
        '+ "limit time range"=0 "internal seed"=0 "include delays"=0\n'
        '+ "frequency at max power"=1 "input signal index"=1\n'
        '+ "stop time"=1 "start time"=0 sensitivity=100f\n'
        '+ "convert noise bins"=1\n'
        )))


    def gen_lsf(self, s):
        self.lsffile.write('{}\n'.format(s))
        return s.split()[0]

    # TODO this should be recursive
    def component_codegen(self, c, parent=None, model=None):
        # In interconnect, we need to catch some subclasses first to
        # make sure they are not generated as they are parent classes.
        # This is a subtle but still ugly workaround for now.
        if isinstance(c, rm.Ammeter):
            tmp_name = self.add_osc(c, parent)
        elif isinstance(c, rm.ConnectionPoint):
            # This is what it is going to look like for Electrical
            # for item in c.components():
                # self.component_codegen(item)
            # tmp_name = c.nodename
            # tmp_name = self.add_sparam(c)
            sch_x, sch_y = self.node_sch_coord(c.coord)
            # note: positions and sch_* have 200 times difference

            ooscs = model.nodes[c.coord[0]][c.coord[1]].ammeters


            conns = ''

            dir_to_port = { 'E':'1', 'W':'2', 'N':'3', 'S':'4' }
            dir_to_port_fiber = { 'E':'1', 'W':'2', 'N':'2', 'S':'1' }

            for d, oosc in ooscs.items():
                conns += self.conn_frmt.format(
                                    i=c.uid,
                                    this_port='port '+dir_to_port[d],
                                    other=oosc.sname,
                                    other_port='input')

            for d,p in dir_to_port.items():
                l = model.get_adjacent_link(parent, d=d)
                if l is not None:
                    conns += self.conn_frmt.format(
                            i=c.uid,
                            this_port='port '+dir_to_port[d],
                            other=l.resistance.sname,
                            other_port='port '+dir_to_port_fiber[d])


            ret = self.gen_lsf(self.sparam_lsf_format.format(i=c.uid,
                                                  sch_x=200*sch_x,
                                                  sch_y=200*sch_y,
                                                  conns=conns))
                                                  
            tmp_name=ret
        elif isinstance(c, rm.VoltageSource):
            # tmp_name = self.add_v2(c)
            tmp_name = self.add_source_compound(c, model)
        elif isinstance(c, rm.Resistance):
            tmp_name = self.add_r2(c, parent)
        else:
            tmp_name = ''
            print("error")
        c.sname = tmp_name

    def add_source_compound(self, c, model):

        size = len(model.mesh)
        pos = self.nodename_to_coord(c.node1)
        sch_x, sch_y = self.node_sch_coord(pos)
        
        if sch_x == 0:
            sch_x -= sch_offset*3
        elif sch_x == size-1:
            sch_x += sch_offset*3
        elif sch_y == 0:
            sch_y -= sch_offset*8
        elif sch_y == size-1:
            sch_y += sch_offset*8

        # the ring resonator at the entry point
        ret = self.gen(self.__ringfrmt.format(i=self.counters['ring'],
                                              nn1=c.node1+'_04',
                                              nn2=c.node1,
                                              sch_x=sch_x,
                                              sch_y=-sch_y,
                                              custom=''))
        self.counters['ring'] += 1

        # first Y branch
        self.gen(self.__yfrmt.format(i=self.counters['y'],
                                     nn1=c.node1+'_04',
                                     nn2=c.node1+'_05',
                                     nn3=c.node1+'_03',
                                     sch_x=sch_x-1,
                                     sch_y=-sch_y,
                                     custom='sch_r=180'))
        self.counters['y'] += 1

        # second Y branch
        self.gen(self.__yfrmt.format(i=self.counters['y'],
                                     nn1=c.node1+'_03',
                                     nn2=c.node1+'_02',
                                     nn3=c.node1+'_01',
                                     sch_x=sch_x-2,
                                     sch_y=-sch_y+1,
                                     custom='sch_r=180'))
        self.counters['y'] += 1

        self.gen(self.__v2frmt.format(i=self.counters['cwl'],
                                      nn1=c.node1+'_01',
                                      sch_x=sch_x-3,
                                      sch_y=-sch_y+2,
                                      custom=''))
        self.counters['cwl'] += 1

        self.gen(self.__pwmfrmt.format(i=self.counters['pwm'],
                                       nn1=c.node1+'_02',
                                       sch_x=sch_x-3,
                                       sch_y=-sch_y,
                                       custom=''))
        self.counters['pwm'] += 1

        self.gen(self.__r2frmt.format(i=self.counters['fiber'],
                                      r=0.1,
                                      uname='',
                                      nn1=c.node1+'_06',
                                      nn2=c.node1+'_05',
                                      sch_x=sch_x-2,
                                      sch_y=-sch_y-1,
                                      custom=''))
        self.counters['fiber'] += 1

        self.gen(self.__yfrmt.format(i=self.counters['y'],
                                     nn1=c.node1+'_07',
                                     nn2=c.node1+'_06',
                                     nn3=c.node1+'_08',
                                     sch_x=sch_x-3,
                                     sch_y=-sch_y-2,
                                     custom=''))
        self.counters['y'] += 1

        self.gen(self.__pwmfrmt.format(i=self.counters['pwm'],
                                       nn1=c.node1+'_07',
                                       sch_x=sch_x-4,
                                       sch_y=-sch_y-2,
                                       custom=''))
        self.counters['pwm'] += 1

        self.gen(self.__yfrmt.format(i=self.counters['y'],
                                     nn1=c.node1+'_08',
                                     nn2=c.node1+'_09',
                                     nn3=c.node1+'_10',
                                     sch_x=sch_x-2,
                                     sch_y=-sch_y-3,
                                     custom=''))
        self.counters['y'] += 1

        self.gen(self.__pwmfrmt.format(i=self.counters['pwm'],
                                       nn1=c.node1+'_10',
                                       sch_x=sch_x-1,
                                       sch_y=-sch_y-4,
                                       custom=''))
        self.counters['pwm'] += 1

        self.gen(self.__v2frmt.format(i=self.counters['cwl'],
                                      nn1=c.node1+'_09',
                                      sch_x=sch_x-1,
                                      sch_y=-sch_y-2,
                                      custom=''))
        self.counters['cwl'] += 1

        return ret

        
    def ic_codegen(self, node, val):
        self.gen(self.__icfrmt.format(node=node, val=val))

    def gen(self, s):
        if not self.cache_only:
            self.file.write('{}\n'.format(s))
        return s.split()[0]



