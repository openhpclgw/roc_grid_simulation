#! /usr/bin/env python
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import *


# in COMSOL I still can't create a 2D region that is to be taken as
# fixed temperature in our case. Which may not make sense with
# physical/realistic use cases, anyway.
# So, here we are creating a heat problem that includes source and sink
# regions, but then when comparing against comsol data, we need to get
# the center slice which doesn't have any source or sink within.
# so we need to create a heat problem larger than the comsol data, so
# that we can get the center slice that has the same data resolution
def main():
    size = int(sys.argv[1])
    mesh_size = int(sys.argv[2])

    comp_size = int(size/(mesh_size-2))
    if comp_size != size/(mesh_size-2):
        print('Divisibility error')

    prb_size = size+2*comp_size  # one component on each side

    source = (0, 0, comp_size, prb_size)
    sink = [(comp_size, 0, prb_size-comp_size, comp_size),
            (prb_size-comp_size, 0, comp_size, prb_size),
            (comp_size, prb_size-comp_size, prb_size-comp_size, comp_size)]


    conductance = 10**-3  # this'll be used as resistance directly
    hp = HeatProblem(prb_size, source, sink, conductance, src_val=1.)

    m = ROCModel(mesh_size)
    m.load_problem(hp)
    m.run_spice_solver(filename='problem')

    from_mesh = m.final_grid[comp_size:comp_size+size,
                             comp_size:comp_size+size]
    from_comsol = load_grid_from_comsol_csv('comsol_data/left_src_others_sink_comsol.csv')

    img_name = '{}_'+str(size)+'_m'+str(mesh_size)
    
    # plot_heatmap(m, filename=img_name.format('heatmap'))
    # plot_surface(m, filename=img_name.format('surface'))
    plot_errmap(from_comsol, from_mesh, filename=img_name.format('err'))



    # rect = 0.1,0.2,0.8,0.7
    # fig = plt.figure(figsize=(15,5))
    # ax = fig.add_axes(rect)

    # img_name = 'cross_comp_s'+str(size)+'_m'+str(mesh_size)+'_dir{}'

    # init_ax(ax, x_lim=size, y_lim=1)
    # mesh_sect = get_isect(from_mesh, 150)
    # comsol_sect = get_isect(from_comsol, 150)
    # custom_plot(ax, range(0,300), mesh_sect, label='Electrical Mesh')
    # custom_plot(ax, range(0,300), comsol_sect, label='Comsol')

    # # ugh
    # plt.legend(handles=datasets)

    # do_plots(img_name.format('i'))

    # ax.clear()

    # init_ax(ax, x_lim=size, y_lim=1)
    # mesh_sect = get_jsect(from_mesh, 150)
    # comsol_sect = get_jsect(from_comsol, 150)
    # custom_plot(ax, range(0,300), mesh_sect, label='Electrical Mesh')
    # custom_plot(ax, range(0,300), comsol_sect, label='Comsol')

    # #ugh
    # plt.legend()

    # do_plots(img_name.format('j'))



# end of main

def get_isect(data, idx):
    i_len = data.shape[0]
    return [r[0] for r in data[0:i_len, idx:idx+1]]

def get_jsect(data, idx):
    j_len = data.shape[1]
    return data[idx:idx+1, 0:j_len][0]

datasets = []
def custom_plot(ax, x, y, label=''):
    ax.set_ylim((0, max(y)*1.1))
    handle, = ax.plot(x, y, label=label,
                            marker='o',
                            markerfacecolor='none',
                            markeredgewidth=2)

    datasets.append(handle)

def init_ax(ax, x_lim, y_lim):
    num_x_ticks=10
    # put only 8 ticks on x axis
    ticks = [i for i in range(0, x_lim, int(x_lim/num_x_ticks))]

    ax.set_xlabel('Grid point')
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)

    ax.set_ylabel('Potential/Temperature')

    ax.set_xlim(0, x_lim)

    ax.grid(b=True, axis='x')
    ax.grid(b=True, axis='y', linestyle='dashed')

    ax.spines['top'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)

def do_plots(filename=''):
    if filename != '':
        plt.savefig(filename)
        plt.savefig(filename+'.eps', format='eps', dpi=1000)
    else:
        plt.show()
    

if __name__=='__main__':
    main()
