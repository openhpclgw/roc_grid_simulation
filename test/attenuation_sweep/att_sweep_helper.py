import sys
from roc_model import ROCModel
import matplotlib.pyplot as plt
from heat_problem import HeatProblem
from analysis_utils import *
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('mesh_size', type=int)
parser.add_argument('multiplier', type=int, default=10000)
parser.add_argument('--generate', action='store_true')
parser.add_argument('--include-boundary-conditions', action='store_true')
parser.add_argument('--save-to-file', action='store_true')

args = parser.parse_args()

mesh_size = args.mesh_size
multiplier = args.multiplier
save_to_file = args.save_to_file
include_boundary_conditions = args.include_boundary_conditions
gen_interconnect_script = args.generate
base_att = int(0.001*multiplier)

def att_sweep(hp, scr_name, out_filename, working_dir):
    
    # quick workaround. size refers to problem size
    size=mesh_size

    m1 = ROCModel(mesh_size)
    m1.load_problem(hp)
    m1.run_interconnect_solver(filename=scr_name,
                              gen_script=gen_interconnect_script,
                              get_results=not gen_interconnect_script)

    if not gen_interconnect_script:
        m2 = ROCModel(mesh_size, sidelinks=False)
        m2.load_problem(hp)
        m2.run_spice_solver()

        normed_m1_grid = normalize_grid(m1.final_grid)
        normed_m2_grid = normalize_grid(m2.final_grid)
        num_grid, num_steps = hp.numerical_solve(num_steps=10000)
        optonum_grid, optonum_steps = hp.numerical_solve(num_steps=10000,
                                          initial_values=normed_m1_grid)

        plot_heatmap_from_grid(normed_m1_grid,
                               filename=out_filename.format('opt'))
        plot_heatmap_from_grid(normed_m2_grid,
                               filename=out_filename.format('elec'))
        plot_heatmap_from_grid(optonum_grid,
                               filename=out_filename.format('optonum'))
        plot_heatmap_from_grid(num_grid,
                               filename=out_filename.format('num'))
        plot_errmap(normed_m1_grid, normed_m2_grid, lim=1.0,
                    filename=out_filename.format('opt_minus_elec'))
        plot_errmap(optonum_grid, normed_m2_grid, lim=1.0,
                    filename=out_filename.format('optonum_minus_elec'))
        plot_errmap(num_grid, normed_m2_grid, lim=1.0,
                    filename=out_filename.format('num_minus_elec'))

        err = np.absolute(normed_m1_grid-normed_m2_grid)
        optonum_err = np.absolute(optonum_grid-normed_m2_grid)
        num_err = np.absolute(num_grid-normed_m2_grid)
 
        max_err = np.absolute(normed_m1_grid-normed_m2_grid).max()
        optonum_max_err = np.absolute(optonum_grid-normed_m2_grid).max()
        num_max_err = np.absolute(num_grid-normed_m2_grid).max()

        rect = 0.1,0.2,0.8,0.7
        fig = plt.figure(figsize=(15,5))
        ax = fig.add_axes(rect)

        init_ax(ax, x_lim=size, y_lim=1)
        custom_plot(ax, range(0,size), normed_m1_grid[int(size/2),:],
                    label='Optical')
        custom_plot(ax, range(0,size), normed_m2_grid[int(size/2),:],
                    label='Electrical')
        custom_plot(ax, range(0,size), optonum_grid[int(size/2),:],
                    label='Optic+numeric')
        custom_plot(ax, range(0,size), num_grid[int(size/2),:],
                    label='Numeric')

        plt.legend(handles=datasets)
        do_plots(filename=out_filename.format('opt_vs_elec_midrow'))

        ax.clear()
        
        data1 = normed_m1_grid[:,int(size/2)]
        data2 = normed_m2_grid[:,int(size/2)]
        # y_lim = max(data1.max(), data2.max())
        y_lim = 1.0
        print(y_lim)
        init_ax(ax, x_lim=size, y_lim=y_lim)
        custom_plot(ax, range(0,size), normed_m1_grid[:,int(size/2)],
                    label='Optical')
        custom_plot(ax, range(0,size), normed_m2_grid[:,int(size/2)],
                    label='Electrical')
        custom_plot(ax, range(0,size), optonum_grid[:,int(size/2)],
                    label='Optic+numeric')
        custom_plot(ax, range(0,size), num_grid[:,int(size/2)],
                    label='Numeric')

        plt.legend()
        do_plots(filename=out_filename.format('opt_vs_elec_midcol'))

        except_idxs = m2.src_idxs | m2.snk_idxs

        print('Report')
        if save_to_file:
            with open(working_dir+'errors', 'a') as f:
                f.write('{} {} {}\n'.format(max_err,
                               mean_except_idxs(err, except_idxs),
                               median_except_idxs(err, except_idxs)))
        else:
            print('Max error : ', max_err)
            print('Mean error : ', mean_except_idxs(err, except_idxs))
            print('Median error : ',
                    median_except_idxs(err, except_idxs))

            print('Max optonum error : ', optonum_max_err)
            print('Mean optonum error : ', 
                        mean_except_idxs(optonum_err, except_idxs))
            print('Median optonum error : ',
                    median_except_idxs(optonum_err, except_idxs))
            print('Optonum steps : ', optonum_steps)

            print('Max num error : ', num_max_err)
            print('Mean num error : ', 
                        mean_except_idxs(num_err, except_idxs))
            print('Median num error : ',
                    median_except_idxs(num_err, except_idxs))
            print('Num steps : ', num_steps)

def mean_except_idxs(data, except_idxs):
    if include_boundary_conditions:
        return np.mean(data)
    else:
        if not isinstance(except_idxs, set):
            print('Except_idxs must be a set')

        val_sum = 0.
        val_count = 0

        for idx, val in np.ndenumerate(data):
            if idx not in except_idxs:
                val_sum += val
                val_count += 1

        return val_sum/val_count

def median_except_idxs(data, except_idxs):
    if include_boundary_conditions:
        return np.median(data)
    else:
        if not isinstance(except_idxs, set):
            print('Except_idxs must be a set')

        valid_values = []
        for idx, val in np.ndenumerate(data):
            if idx not in except_idxs:
                valid_values.append(val)

        valid_values.sort()
        val_count = len(valid_values)

        if val_count % 2 == 1:
            return valid_values[val_count//2+1]
        else:
            return (valid_values[val_count//2]+
                    valid_values[val_count//2+1])/2


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

