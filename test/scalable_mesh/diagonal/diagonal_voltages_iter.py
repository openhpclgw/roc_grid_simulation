#! /usr/bin/env python
import sys
import numpy as np
import itertools as it
import matplotlib as mpl
import matplotlib.pyplot as plt
from roc_model import ROCModel
from heat_problem import HeatProblem
from analysis_utils import (aggregate_current_vectors,
                            print_current_table,
                            energy_flow,
                            plot_errmap,
                            plot_heatmap_from_grid,
                            plot_heatmap)

filename='tmp/diag_v_{}'

max_exp_prob_size = 7
max_prob_size = 2**max_exp_prob_size+1
def main():

    use_cached = True 

    rect = 0.1,0.2,0.8,0.7
    fig = plt.figure(figsize=(15,5))
    ax = fig.add_axes(rect)

    diag_data = {}
    full_data = {}
    interp_data = {}
    scale_data = {}

    for exp_size in range(2, max_exp_prob_size+1):
        prob_size = 2**exp_size+1
        comp_size = 1
        source = (0, 0, comp_size, comp_size)
        sink = (prob_size-comp_size, 
                prob_size-comp_size,
                comp_size, comp_size)
        cond_exp = -3
        conductance = 10**cond_exp
        hp = HeatProblem(prob_size, source, sink, conductance, src_val=1.)
        mesh = ROCModel(prob_size)

        res_grid = get_results(mesh, hp, cached=use_cached)
        full_data[prob_size] = res_grid

        prob_range = range(prob_size)
        tmp_data = [res_grid[i,i] for i in prob_range]
        diag_data[prob_size] = tmp_data

        interp_data[prob_size] = get_interp2d(res_grid, max_prob_size)
        scale_data[prob_size] = get_scale2d(res_grid, max_prob_size)

    # we have the raw data in the dictionary now
    # let's do some plots

    print('Finished reading data')

    # main plot
    for exp_size in range(2, max_exp_prob_size+1):
        prob_size = 2**exp_size+1
        exp_prob_range = range(0, max_prob_size,
                                  int(2**(max_exp_prob_size-exp_size)))

        custom_plot(ax, exp_prob_range, diag_data[prob_size],
                    label='{}-by-{}'.format(prob_size, prob_size))

    plt.legend()
    do_plot('scale_comparison')
    ax.clear()

    print('Finished main plot')

    # point errors
    base_data = diag_data[max_prob_size]
    for exp_size in range(2, max_exp_prob_size):
        prob_size = 2**exp_size+1
        exp_prob_range = range(0, max_prob_size,
                                  int(2**(max_exp_prob_size-exp_size)))

        interp_res_grid = get_interp1d(diag_data[prob_size], max_prob_size)

        prob_range = range(max_prob_size)
        custom_plot(ax, prob_range, interp_res_grid-base_data,
                    label='{}-by-{}'.format(prob_size, prob_size))



    plt.legend()
    do_plot('error_vs_largest')
    ax.clear()

    print('Finished point errors')

    # average absolute error
    base_data = full_data[max_prob_size]
    aae_interp = {}
    aae_scale = {}
    aae_pwise = {}
    for exp_size in range(2, max_exp_prob_size):
        prob_size = 2**exp_size+1
        # exp_prob_range = range(0, max_prob_size,
                                  # int(2**(max_exp_prob_size-exp_size)))

        # interp_res_grid = get_interp2d(full_data[prob_size],
                                       # max_prob_size)
        sum_err = abs(interp_data[prob_size]-base_data).sum()
        aae_interp[prob_size] = sum_err/(max_prob_size**2)

        # scale_res_grid = get_scale2d(full_data[prob_size],
                                     # max_prob_size)
        sum_err = abs(scale_data[prob_size]-base_data).sum()
        aae_scale[prob_size] = sum_err/(max_prob_size**2)

        factor = int((max_prob_size-1)/(prob_size-1))
        my_fdata = full_data[prob_size]
        raw_data = [abs(my_fdata[i][j]-base_data[i*factor][j*factor]) 
                      for i,j in it.product(range(prob_size),
                                            range(prob_size))]
        sum_err = sum(raw_data)
        aae_pwise[prob_size] = sum_err/(prob_size**2)


    custom_plot(ax, [k for k,v in aae_pwise.items()],
                    [v for k,v in aae_pwise.items()],
                ticks=[2**s+1 for s in range(2, max_exp_prob_size+1)],
                label='Pointwise')
    custom_plot(ax, [k for k,v in aae_interp.items()],
                    [v for k,v in aae_interp.items()],
                ticks=[2**s+1 for s in range(2, max_exp_prob_size+1)],
                label='Interpolated')
    custom_plot(ax, [k for k,v in aae_scale.items()],
                    [v for k,v in aae_scale.items()],
                ticks=[2**s+1 for s in range(2, max_exp_prob_size+1)],
                label='Scaled')

    plt.legend()
    do_plot('avg_abs_err')

    print('Finished absolute errors')

    # error maps
    for exp_size in range(2, max_exp_prob_size):
        prob_size = 2**exp_size+1

        interp_res_grid = interp_data[prob_size]
        plot_heatmap_from_grid(interp_res_grid,
                               filename='heatmap_interp_{}'.format(
                                                             prob_size))
        plot_errmap(interp_res_grid, base_data,
                    filename='errmap_interp{}'.format(prob_size))

        scale_res_grid = scale_data[prob_size]
        plot_heatmap_from_grid(scale_res_grid,
                               filename='heatmap_scale_{}'.format(
                                                            prob_size))
        plot_errmap(scale_res_grid, base_data,
                    filename='errmap_scale{}'.format(prob_size))

    print('Finished heat and error maps')

def get_results(mesh, hp, cached=False):
    import os.path

    mesh.load_problem(hp)
    f = filename.format(mesh.w)
    exists = os.path.exists(f+'.out')
    if cached and exists:
        print('\t\tUsing cache ' + f)
        mesh.init_from_cache(f)
    else:
        mesh.run_spice_solver(f)
    return mesh.final_grid

def custom_plot(ax, x, y, label='', ticks=None):
    ax.plot(x, y, label=label,
                  marker='o',
                  markerfacecolor='none',
                  markeredgewidth=2)

    if ticks is None:
        prob_size = len(x)
        # put only 8 ticks on x axis
        ticks = [i for i in range(0, prob_size, int(prob_size/8) if
            prob_size >= 8 else 1)]

    else:
        prob_size = max([int(t) for t in ticks])

    ax.set_xlabel('Grid point')
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)

    ax.set_ylabel('Potential')

    # ax.set_ylim((min(y), max(y)*1.1))
    ax.set_xlim(0,prob_size-1)

    ax.grid(b=True, axis='x')
    ax.grid(b=True, axis='y', linestyle='dashed')

    ax.spines['top'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)

def get_interp1d(grid, interp_size):

    def get_interpolator():
        import scipy.interpolate
        grid_size = len(grid)

        assert interp_size >= grid_size
        
        _tmp = (interp_size-1)/(grid_size-1)
        interp_factor = _tmp

        interp_vals = [i*interp_factor for i in range(0, grid_size)]

        return scipy.interpolate.interp1d(interp_vals,
                                          grid)

    interpolator = get_interpolator()
    # res = np.zeros(interp_size)
    # for i in range(interp_size):
        # res[i] = interpolator(i)

    # return res

    return interpolator(range(interp_size))

def get_interp2d(grid, interp_size):

    def get_interpolator():
        import scipy.interpolate
        assert grid.shape[0] == grid.shape[1]
        grid_size = grid.shape[0]

        assert interp_size >= grid_size
        
        _tmp = (interp_size-1)/(grid_size-1)
        interp_factor = _tmp

        interp_vals = [i*interp_factor for i in range(0, grid_size)]
        return scipy.interpolate.interp2d(interp_vals,
                                          interp_vals,
                                          grid)

    interpolator = get_interpolator()
    res = np.zeros((interp_size, interp_size))

    for i,j in it.product(range(interp_size), range(interp_size)):
        res[i,j] = interpolator(i,j)

    return res

def get_scale2d(grid, scale_size):
    assert grid.shape[0] == grid.shape[1]
    grid_size = grid.shape[0]
    assert scale_size >= grid_size
    
    scale_factor = scale_size/grid_size
    res = np.zeros((scale_size, scale_size))

    for i,j in it.product(range(scale_size), range(scale_size)):
        res[i,j] = grid[int(i/scale_factor), int(j/scale_factor)]

    return res

def do_plot(filename=''):
    if filename == '':
        plt.show()
    else:
        plt.savefig(filename+'.png')
        plt.savefig(filename+'.eps')

if __name__=='__main__':
    main()
