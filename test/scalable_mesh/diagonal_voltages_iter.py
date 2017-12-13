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
                            plot_heatmap)

filename='tmp/diag_v_{}'

max_exp_prob_size = 5
max_prob_size = 2**max_exp_prob_size+1
def main():

    use_cached = True 

    rect = 0.1,0.2,0.8,0.7
    fig = plt.figure(figsize=(15,5))
    ax = fig.add_axes(rect)

    data = {}

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

        interp_res_grid = get_interp1d(
                [res_grid[i,i] for i in range(len(res_grid))],
                max_prob_size)

        prob_range = range(max_prob_size)
        tmp_data = interp_res_grid
        print(tmp_data)

        data[prob_size] = tmp_data
        custom_plot(ax, prob_range, tmp_data,
                    label='{}-by-{}'.format(prob_size, prob_size))

        prob_range = range(prob_size)
        tmp_data = [res_grid[i,i] for i in prob_range]
        data[prob_size] = tmp_data
        exp_prob_range = range(0, max_prob_size,
                                  int(2**(max_exp_prob_size-exp_size)))
        custom_plot(ax, exp_prob_range, tmp_data,
                    label='{}-by-{}'.format(prob_size, prob_size))

    plt.legend()
    plt.show()

    # plt.cla()
    # # ax = fig.add_axes(rect)

    # base = data[max_prob_size]
    # for exp_size in range(2, max_exp_prob_size):
        # prob_size = 2**exp_size
        # cur_data = data[prob_size]
        # factor = int(max_prob_size/prob_size)
        # err = [cur_data[i]-base[i*factor] for i in range(len(cur_data))]

        # exp_prob_range = range(0, 2**max_exp_prob_size,
                                  # int(2**max_exp_prob_size/prob_size))
        # print(err)
        # custom_plot(ax, [i for i in exp_prob_range], err,
                    # label='{}-by-{}'.format(prob_size, prob_size))

    # plt.legend()
    # plt.show()

def get_results(mesh, hp, cached=False):
    import os.path

    mesh.load_problem(hp)
    f = filename.format(mesh.w)
    exists = os.path.exists(f+'.out')
    if cached and exists:
        print('Using cache' + f)
        mesh.init_from_cache(f)
    else:
        mesh.run_spice_solver(f)
    return mesh.final_grid

def custom_plot(ax, x, y, label=''):
    ax.plot(x, y, label=label,
                  marker='o',
                  markerfacecolor='none',
                  markeredgewidth=2)

    prob_size = len(x)
    # put only 8 ticks on x axis
    ticks = [i for i in range(0, prob_size, int(prob_size/8) if
        prob_size >= 8 else 1)]

    ax.set_xlabel('Grid point')
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)

    ax.set_ylabel('Potential')

    ax.set_ylim((min(y), max(y)*1.1))
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
        # interp_range = range(0, grid_size, interp_factor)

        interp_vals = [i*interp_factor for i in range(0, grid_size)]
        # interp_vals = np.linspace(0, interp_size, interp_factor)
        print(interp_vals)

        return scipy.interpolate.interp2d(interp_vals,
                                          interp_vals,
                                          grid)

    interpolator = get_interpolator()
    res = np.zeros((interp_size, interp_size))

    for i,j in it.product(range(interp_size), range(interp_size)):
        res[i,j] = interpolator(i,j)

    return res


if __name__=='__main__':
    main()