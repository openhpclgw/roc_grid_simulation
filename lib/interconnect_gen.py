import re
import roc_model as rm
import itertools as it
import numpy as np
import subprocess as sp
from common import *

coord_to_sch_ratio = 20.
mid_distance = coord_to_sch_ratio/2
q_distance = mid_distance/2
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
        self.__v2frmt = 'X_CWL_'+cg+' {nn1} \"CW Laser\" \"internal seed\"=7434 power={power} '+self.__schfrmt
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
                                  'set("s parameters filename","{f}.txt");\n'
                                  'set("x position", {sch_x:4.2f});\n'
                                  'set("y position", {sch_y:4.2f});\n'
                                  '{conns}\n')
        self.conn_frmt = ('connect("{i}","{this_port}",'
                                  ' "{other}", "{other_port}");\n')
        self.__pwmfrmt = 'X_OPWM_'+cg+' {nn1} \"Optical Power Meter\"'+self.__schfrmt
        self.__yfrmt = 'X_Y_'+cg+' {nn1} {nn2} {nn3} \"Waveguide Y Branch\"'+self.__schfrmt
        self.__ringfrmt = 'X_RING_'+cg+' {nn1} {nn2} \"Single Bus Ring Resonator\" frequency=230T '+self.__schfrmt
        self.ring_lsf_format = ('addelement("Single Bus Ring Resonator");\n'
                                'set("name","RING'+cg+'{ringdir}");\n'
                                'set("x position", {sch_x:4.2f});\n'
                                'set("y position", {sch_y:4.2f});\n'
                                '{conns}\n')

    def create_script(self, roc_model, suffix=''):
        if not self.cache_only:
            self.file = open(self.rel_in_path(suffix), 'w')
            self.lsffile = open(self.rel_lsf_path(suffix), 'w')
        # mesh size
        if roc_model.norton:
            self.mesh_size = len(roc_model.mesh)+1
        else:
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
        # print(full_range)
        self.add_block_comment("Node subcircuits")
        for i, j in it.product(full_range, full_range):
            n = roc_model.nodes[i][j]
            self.add_comment("Node " + str((i, j)))
            for c in n.components():
                self.component_codegen(c, n, model=roc_model)

            # generate the initial condition
            if n.ic != 0.:
                self.ic_codegen(n.sname, n.ic)

        if not roc_model.norton:
            # generate source
            self.add_block_comment("Sources")
            for v in roc_model.src:
                self.component_codegen(v, model=roc_model)

            # generate sink
            self.add_block_comment("Sinks")
            for s in roc_model.snk:
                self.component_codegen(s, model=roc_model)
        else:
            # generate loop components for Norton circuit
            self.add_block_comment('Loop Components')
            for loop in roc_model.iter_loops():
                self.add_comment('Loop ' + str(loop.coord))
                for mr in loop.mesh_resistances():
                    self.add_comment('Loop resistance '+mr.node1+' '+
                                                        mr.node2)
                    for c in mr.components():
                        # print('loop res ', mr.node1, mr.node2)
                        self.component_codegen(c, model=roc_model, parent=mr)
                        # print('Sname = ', c.sname)

                        if isinstance(c, rm.Resistance):
                            left_coord = self.nodename_to_coord(mr.node1)
                            right_coord = self.nodename_to_coord(mr.node2)

                            print('Mesh resistance btw ',
                                  left_coord, right_coord)

                            left_dir = mr.node1[-1]
                            right_dir = mr.node2[-1]

                            left_conn_id = self.cid_gr.format(i=roc_model.nodes[left_coord[0]][left_coord[1]].conn_point.uid)
                            right_conn_id = self.cid_gr.format(i=roc_model.nodes[right_coord[0]][right_coord[1]].conn_point.uid)

                            
                            self.gen_lsf(self.conn_frmt.format(
                                            i=c.sname[0],
                                            this_port='port 1',
                                            other='RING'+left_conn_id+left_dir,
                                            other_port='port 2'))

                            self.gen_lsf(self.conn_frmt.format(
                                            i=c.sname[1],
                                            this_port='port 1',
                                            other='RING'+right_conn_id+right_dir,
                                            other_port='port 2'))
                for c in loop.other_components():
                    self.add_comment('Loop component')
                    self.component_codegen(c, model=roc_model, parent=loop)

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
        tmp_name = name
        if name[-1:].isalpha():
            tmp_name = name[:-1]
        coord_str = tmp_name[1:].split('_')
        print('Nodename = ', name, ' Coord = ',
              str((int(coord_str[0]), int(coord_str[1]))))
        return Coord(int(coord_str[0]), int(coord_str[1]))

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

        if roc_model.norton:
            for loop in roc_model.iter_loops():
                for c in loop.curmeters():
                    c.current = result_dict[c.sname]

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

    def add_r2(self, r, model=None, parent=None):
        custom = ''
        sch_x, sch_y = 0, 0

        if parent is not None:
            if not isinstance(parent, rm.NortonLoop):
                node1_sch = self.node_sch_coord(parent.nodeblock1.coord)
                node2_sch = self.node_sch_coord(parent.nodeblock2.coord)
                if parent.orientation == 'V':
                    custom = '"rotated"=true'
            else:
                node1_sch = self.node_sch_coord(parent.coord)
                node2_sch = self.node_sch_coord(parent.coord)

            sch_x, sch_y = self.midpoint(node1_sch, node2_sch)

        ret = self.create_link_compound('R'+str(r.uid), sch_x, sch_y,
                                        attenuation=r.r,
                                        power_left=0.,
                                        power_right=0.)

        if isinstance(parent, rm.NortonLoop):
            assert False
            pos = self.nodename_to_coord(r.node2)
            if pos[1] == 0:
                ringdir = 'W'
            else:
                ringdir = 'E'

            self.gen_lsf(self.conn_frmt.format(i=ret[1],
                this_port='port 1',
                other='RING'+
                    self.cid_gr.format(i=model.nodes[pos[0]][pos[1]].conn_point.uid)+ ringdir,
                other_port='port 2'))

        self.r_counter += 1
        return ret

    def conn_uid_str(self, model, coord):
        return self.cid_gr.format(i=model.nodes[coord[0]][coord[1]].conn_point.uid)

    def add_current_source(self, cs, model=None, parent=None):

        assert model is not None
        assert parent is not None

        custom = ''
        sch_x, sch_y = 0, 0

        n1 = self.nodename_to_coord(cs.node1)
        n2 = self.nodename_to_coord(cs.node2)
        pos = parent.coord

        node1_sch = self.node_sch_coord(n1)
        node2_sch = self.node_sch_coord(n2)

        sch_x, sch_y = self.midpoint(node1_sch, node2_sch)

        right_power = False

        if n1.is_h(n2) and n1.i != 0: # bottom
            right_power = True
        if (not n1.is_h(n2)) and n1.j != 0: # right
            right_power = True

        if right_power:
            pl, pr = 0, cs.i*4
        else:
            pr, pl = 0, cs.i*4


        ret = self.create_link_compound('Cur_Source'+str(cs.uid),
                                        sch_x, sch_y,
                                        attenuation=0,
                                        power_left=pl,
                                        power_right=pr)

        self.r_counter += 1


        # print('Node coords ', n1, n2)

        if n1.is_h(n2):
            if n1.i == 0:
                # source is on top of the loop (and the whole mesh)
                self.gen_lsf(self.conn_frmt.format(i=ret[0],
                  this_port='port 1',
                  other='RING'+self.conn_uid_str(model,(pos[0],pos[1]))+'E',
                  other_port='port 2'))

                self.gen_lsf(self.conn_frmt.format(i=ret[1],
                  this_port='port 1',
                  other='RING'+self.conn_uid_str(model,(pos[0],pos[1]+1))+'W',
                  other_port='port 2'))

            else:
                # it must be on the bottom
                self.gen_lsf(self.conn_frmt.format(i=ret[0],
                  this_port='port 1',
                  other='RING'+self.conn_uid_str(model,(pos[0]+1,pos[1]))+'E',
                  other_port='port 2'))

                self.gen_lsf(self.conn_frmt.format(i=ret[1],
                  this_port='port 1',
                  other='RING'+self.conn_uid_str(model,(pos[0]+1,pos[1]+1))+'W',
                  other_port='port 2'))
        else:
            if n1.j == 0:
                # source is on left of the loop (and the whole mesh)
                self.gen_lsf(self.conn_frmt.format(i=ret[0],
                  this_port='port 1',
                  other='RING'+self.conn_uid_str(model,(pos[0]+1,pos[1]))+'N',
                  other_port='port 2'))

                self.gen_lsf(self.conn_frmt.format(i=ret[1],
                  this_port='port 1',
                  other='RING'+self.conn_uid_str(model,(pos[0],pos[1]))+'S',
                  other_port='port 2'))

            else:
                # it must be on the right
                self.gen_lsf(self.conn_frmt.format(i=ret[0],
                  this_port='port 1',
                  other='RING'+self.conn_uid_str(model,(pos[0]+1,pos[1]+1))+'N',
                  other_port='port 2'))

                self.gen_lsf(self.conn_frmt.format(i=ret[0],
                  this_port='port 1',
                  other='RING'+self.conn_uid_str(model,(pos[0],pos[1]+1))+'S',
                  other_port='port 2'))

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
                                         power=v.v,
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
                    sch_x += sch_offset*1.1

                conn_p = parent.resistance.sname[3]

            elif isinstance(parent, rm.NodeBlock) or isinstance(parent, rm.NortonLoop):
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

                conn_p = osc.node2

        else:
            sch_x, sch_y = 0, 0

        ret = ''
        ret = self.gen(
                self.__oscfrmt.format(i=osc.uid,
                                     uname='',
                                     # nn1=osc.node1,
                                     nn2=conn_p,
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
        '+ frequency=193.1T "reference power"=0.001 RIN=-140\n'
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
        import re
        self.lsffile.write('{}\n'.format(s))

        reobj = re.compile('set\("name","(.*)"\);')
        for line in s.splitlines():
            m = reobj.match(line)
            if m:
                return m.group(1)

    # TODO this should be recursive
    def component_codegen(self, c, parent=None, model=None):
        # In interconnect, we need to catch some subclasses first to
        # make sure they are not generated as they are parent classes.
        # This is a subtle but still ugly workaround for now.
        if isinstance(c, rm.CurrentMeter):
            tmp_name = self.add_osc(c, parent)
        elif isinstance(c, rm.ConnectionPoint): #creates lsf
            # This is what it is going to look like for Electrical
            # for item in c.components():
                # self.component_codegen(item)
            # tmp_name = c.nodename
            # tmp_name = self.add_sparam(c)
            sch_x, sch_y = self.node_sch_coord(c.coord)
            # note: positions and sch_* have 200 times difference

            ooscs = model.nodes[c.coord[0]][c.coord[1]].curmeters


            conns = ''

            dir_to_port = { 'E':'1', 'W':'2', 'N':'3', 'S':'4' }
            dir_to_port_fiber = { 'E':0, 'W':1, 'N':1, 'S':0 }

            print('coord ', c.coord, ' id ', c.uid)

            resonators = {
            'E': self.gen_lsf(self.ring_lsf_format.format(i=str(c.uid),
                                 ringdir='E',
                                 sch_x=200*(sch_x+q_distance),
                                 sch_y=200*sch_y, conns='')),
            'W': self.gen_lsf(self.ring_lsf_format.format(i=str(c.uid),
                                 ringdir='W',
                                 sch_x=200*(sch_x-q_distance),
                                 sch_y=200*sch_y, conns='')),
            'N': self.gen_lsf(self.ring_lsf_format.format(i=str(c.uid),
                                 ringdir='N',
                                 sch_x=200*sch_x,
                                 sch_y=200*(sch_y-q_distance), conns='')),
            'S': self.gen_lsf(self.ring_lsf_format.format(i=str(c.uid),
                                 ringdir='S',
                                 sch_x=200*sch_x,
                                 sch_y=200*(sch_y+q_distance), conns=''))}

            # connection to adjacent links
            for d,p in dir_to_port.items():
                l = model.get_adjacent_link(parent, d=d)
                if l is not None:
                    print('helloo')
                    print('sname ', l.resistance.sname)
                    conns += self.conn_frmt.format(
                            i=resonators[d],
                            this_port='port 2',
                            other=l.resistance.sname[dir_to_port_fiber[d]],
                            other_port='port 1')


            spar = self.gen_lsf(self.sparam_lsf_format.format(i=c.uid,
                                                  sch_x=200*sch_x,
                                                  sch_y=200*sch_y,
                                                  # f=c.get_spar_file(),
                                                  f='spar',
                                                  conns=conns))

            # connection to node oscillators
            for d, oosc in ooscs.items():
                conns += self.conn_frmt.format(
                                    i=spar,
                                    this_port='port '+dir_to_port[d],
                                    other=oosc.sname,
                                    other_port='input')

            for d,p in dir_to_port.items():
                self.gen_lsf(self.conn_frmt.format(i=spar,
                            this_port='port '+str(dir_to_port[d]),
                            other=resonators[d],
                            other_port='port 1'))

            
                                                  
            tmp_name=spar

        elif isinstance(c, rm.BoundaryCond):
            # assert c.v==0
            tmp_name = self.add_v2(c)

        elif isinstance(c, rm.CurrentSource):
            print('Current Source btw ', c.node1, c.node2)
            tmp_name = self.add_current_source(c, model, parent)
        elif isinstance(c, rm.Resistance):
            tmp_name = self.add_r2(c, model, parent)
        else:
            tmp_name = ''
            print(c)
            print("unrecognized type error")
            assert False
        print(tmp_name)
        c.sname = tmp_name

    def add_source_compound(self, c, model):

        size = len(model.mesh)
        pos = self.nodename_to_coord(c.node1)
        print(c.node1, ' > ', pos, ' > ', model.nodes[pos[0]][pos[1]].conn_point.uid)
        sch_x, sch_y = self.node_sch_coord(pos)
        
        if pos[1] == 0:
            sch_x -= mid_distance
        elif pos[1] == size-1:
            sch_x += mid_distance
        elif pos[0] == 0:
            sch_y -= mid_distance
        elif pos[0] == size-1:
            sch_y += mid_distance

        # the ring resonator at the entry point
        # ret = self.gen(self.__ringfrmt.format(i=self.counters['ring'],
                                              # nn1=c.node1+'_04',
                                              # nn2=c.node1,
                                              # sch_x=sch_x,
                                              # sch_y=-sch_y,
                                              # custom=''))
        # self.counters['ring'] += 1

        ret = self.create_link_compound('bc'+str(c.uid), sch_x, sch_y)

        if not model.norton:
            pwm = self.gen(self.__pwmfrmt.format(i=self.counters['pwm'],
                                           nn1=c.node1+'_07',
                                           sch_x=sch_x-4,
                                           sch_y=-sch_y-2,
                                           custom=''))

            self.counters['pwm'] += 1

        # name of the ring resonator must be found
        if not model.norton:
            if pos[1] == 0:
                ringdir = 'W'
            else:
                ringdir = 'E'

            self.gen_lsf(self.conn_frmt.format(i=ret[1],
                this_port='port 1',
                other='RING'+
                    self.cid_gr.format(i=model.nodes[pos[0]][pos[1]].conn_point.uid)+
                    ringdir,
                other_port='port 2'))

            self.gen_lsf(self.conn_frmt.format(i=ret[0],
                                               this_port='port 1',
                                               other=pwm,
                                               other_port='input'))

        else:
            # the source is connected to two nodes in Norton circuits.
            # For now assume that the current source will always be on
            # the west boundary of the loop:

            self.gen_lsf(self.conn_frmt.format(i=ret[0],
                this_port='port 1',
                other='RING'+
                    self.cid_gr.format(i=model.nodes[pos[0]][pos[1]].conn_point.uid)+
                    'S',
                other_port='port 2'))

            self.gen_lsf(self.conn_frmt.format(i=ret[1],
                this_port='port 1',
                other='RING'+
                    self.cid_gr.format(i=model.nodes[pos[0]+1][pos[1]].conn_point.uid)+
                    'S',
                other_port='port 2'))
            


        return ret

    def create_link_compound(self, name, sch_x, sch_y,
                             attenuation=0.,
                             power_left=0.,
                             power_right=0.):


        # first Y branch
        right = self.gen(self.__yfrmt.format(i=self.counters['y'],
                                             nn1=name+'_04',
                                             nn2=name+'_05',
                                             nn3=name+'_03',
                                             sch_x=sch_x-1,
                                             sch_y=-sch_y,
                                             custom='sch_r=180'))
        self.counters['y'] += 1

        # second Y branch
        self.gen(self.__yfrmt.format(i=self.counters['y'],
                                     nn1=name+'_03',
                                     nn2=name+'_02',
                                     nn3=name+'_01',
                                     sch_x=sch_x-2,
                                     sch_y=-sch_y+1,
                                     custom='sch_r=180'))
        self.counters['y'] += 1

        cwl_right = self.gen(self.__v2frmt.format(i=self.counters['cwl'],
                                      nn1=name+'_01',
                                      power=power_right,
                                      sch_x=sch_x-3,
                                      sch_y=-sch_y+2,
                                      custom=''))
        self.counters['cwl'] += 1

        # this measures the current from the right
        self.gen(self.__pwmfrmt.format(i=self.counters['pwm'],
                                       nn1=name+'_02',
                                       sch_x=sch_x-3,
                                       sch_y=-sch_y,
                                       custom=''))
        self.counters['pwm'] += 1

        self.gen(self.__r2frmt.format(i=self.counters['fiber'],
                                      r=attenuation,
                                      uname='',
                                      nn1=name+'_06',
                                      nn2=name+'_05',
                                      sch_x=sch_x-2,
                                      sch_y=-sch_y-1,
                                      custom=''))
        self.counters['fiber'] += 1

        left = self.gen(self.__yfrmt.format(i=self.counters['y'],
                                            nn1=name+'_07',
                                            nn2=name+'_06',
                                            nn3=name+'_08',
                                            sch_x=sch_x-3,
                                            sch_y=-sch_y-2,
                                            custom=''))
        self.counters['y'] += 1

        self.gen(self.__yfrmt.format(i=self.counters['y'],
                                     nn1=name+'_08',
                                     nn2=name+'_09',
                                     nn3=name+'_10',
                                     sch_x=sch_x-2,
                                     sch_y=-sch_y-3,
                                     custom=''))
        self.counters['y'] += 1

        # this measures the current from the left
        self.gen(self.__pwmfrmt.format(i=self.counters['pwm'],
                                       nn1=name+'_10',
                                       sch_x=sch_x-1,
                                       sch_y=-sch_y-4,
                                       custom=''))
        self.counters['pwm'] += 1

        cwl_left = self.gen(self.__v2frmt.format(i=self.counters['cwl'],
                                      nn1=name+'_09',
                                      power=power_left,
                                      sch_x=sch_x-1,
                                      sch_y=-sch_y-2,
                                      custom=''))
        self.counters['cwl'] += 1

        return left, right, name+'_07', name+'_04'

        
    def ic_codegen(self, node, val):
        self.gen(self.__icfrmt.format(node=node, val=val))

    def gen(self, s):
        if not self.cache_only:
            self.file.write('{}\n'.format(s))
        return s.split()[0]



