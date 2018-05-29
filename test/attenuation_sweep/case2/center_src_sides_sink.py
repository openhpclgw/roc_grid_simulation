#! /usr/bin/env python
import sys
from roc_model import ROCModel
import matplotlib.pyplot as plt
from heat_problem import HeatProblem
from analysis_utils import *
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('mesh_size', type=int)
parser.add_argument('multiplier', type=int, default=10000)

args = parser.parse_args()

mesh_size = args.mesh_size
multiplier = args.multiplier
base_att = int(0.001*multiplier)

def main():
    def resistance(size, n1, n2):
        return 0.001*multiplier

    size = int(sys.argv[1])

    N = size
    mesh_size = size
    # source is a 2x2 in the middle
    assert size%2==0
    source = (int(size/2)-1, int(size/2)-1, 2, 2)
    sink = [(0, 0, 1, size),  # left
            (size-1, 0, 1, size),  # right
            (1, 0, size-2, 1),      # top
            (1, size-1, size-2, 1)]   # bottom

    hp = HeatProblem(N, source, sink, resistance)

    working_dir = 'test/attenuation_sweep/case2/size'+str(mesh_size)+'/'
    filename = working_dir+'{}_'+str(size)


    gen_interconnect_script = True
    # assert not gen_interconnect_script


    m1 = ROCModel(mesh_size)
    m1.load_problem(hp)
    m1.run_interconnect_solver(filename=working_dir+'opt_case2_'+str(size)+'_att'+str(base_att),
                              gen_script=gen_interconnect_script,
                              get_results=not gen_interconnect_script)


    if not gen_interconnect_script:
        m2 = ROCModel(mesh_size, sidelinks=False)
        m2.load_problem(hp)
        m2.run_spice_solver()

        normed_m1_grid = normalize_grid(m1.final_grid)
        normed_m2_grid = normalize_grid(m2.final_grid)


        plot_heatmap_from_grid(normed_m1_grid,
                               filename=filename.format('opt'))
        plot_heatmap_from_grid(normed_m2_grid,
                               filename=filename.format('elec'))
        plot_errmap(normed_m1_grid, normed_m2_grid,
                    filename=filename.format('opt_minus_elec'))

        err = np.absolute(normed_m1_grid-normed_m2_grid)
        max_err = np.absolute(normed_m1_grid-normed_m2_grid).max()

        rect = 0.1,0.2,0.8,0.7
        fig = plt.figure(figsize=(15,5))
        ax = fig.add_axes(rect)

        init_ax(ax, x_lim=size, y_lim=1)
        custom_plot(ax, range(0,size), normed_m1_grid[int(size/2),:],
                    label='Optical')
        custom_plot(ax, range(0,size), normed_m2_grid[int(size/2),:],
                    label='Electrical')

        plt.legend(handles=datasets)
        do_plots(filename=filename.format('opt_vs_elec_midrow'))

        ax.clear()

        
        data1 = normed_m1_grid[:,int(size/2)]
        data2 = normed_m2_grid[:,int(size/2)]
        y_lim = max(data1.max(), data2.max())
        print(y_lim)
        init_ax(ax, x_lim=size, y_lim=y_lim)
        custom_plot(ax, range(0,size), normed_m1_grid[:,int(size/2)],
                    label='Optical')
        custom_plot(ax, range(0,size), normed_m2_grid[:,int(size/2)],
                    label='Electrical')


        plt.legend()
        do_plots(filename=filename.format('opt_vs_elec_midcol'))

        print('Report')
        print('Max error : ', max_err)
        print('Mean error : ', np.mean(err))
        print('Median error : ', np.median(err))




# plot functions
datasets = []
def custom_plot(ax, x, y, label=''):
    # ax.set_ylim((0, max(y)*1.1))
    handle, = ax.plot(x, y, label=label,
                            marker='o',
                            markerfacecolor='none',
                            markeredgewidth=2)

    datasets.append(handle)

def init_ax(ax, x_lim, y_lim):
    num_x_ticks=10 if x_lim>=10 else x_lim
    # put only 8 ticks on x axis
    ticks = [i for i in range(0, x_lim, int(x_lim/num_x_ticks))]

    ax.set_xlabel('Grid point')
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)

    ax.set_ylabel('Potential/Temperature')

    ax.set_xlim(0, x_lim-1)
    ax.set_ylim(0, y_lim*1.1)

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
