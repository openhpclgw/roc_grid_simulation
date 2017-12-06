import numpy as np
import itertools as it
import sys
import math


def print_current_table(m, print_checksum=False):
    frs = '{0:>9}    {1:>9}     {2}         {3:>1}'
    print(frs.format('Terminal1', 'Terminal2', 'Current', 'Direction'))
    checksum = 0.
    for mr in m.links:
        abs_current = abs(mr.ammeter.current)
        print(frs.format(str(mr.nodeblock1.coord),
                         str(mr.nodeblock2.coord),
                         '{:>10.6e}'.format(abs_current),
                         mr.cur_direction(mr.ammeter.current)))
        checksum += abs_current

    if print_checksum:
        print('Checksum: {}'.format(checksum))


def print_node_potentials(m):
    frs = '{0:>9}    {1:>12}'
    print(frs.format('Node', 'Potential'))
    for n in m.iter_nodes():
        print(frs.format(str(n.coord), '{: e}'.format(n.potential)))


def print_node_currents(m):
    frs = '{:>9}    {}'
    dict_format = 'E:{E: >14e}   W:{W: >14e}   N:{N: >14e}   S:{S: >14e}'
    print(frs.format('Node', 'Currents (+:in, -:out)'))
    for n in m.iter_nodes():
        print(frs.format(str(n.coord),
                         dict_format.format(**nodal_current_dict(n))))


def energy_flow(m):
    src_out = sum([n.sum_reduce_in_curs() for n in m.snk_nodeblocks()])
    snk_in = sum([n.sum_reduce_out_curs() for n in m.src_nodeblocks()])

    src_in = sum([n.sum_reduce_in_curs() for n in m.src_nodeblocks()])
    snk_out = sum([n.sum_reduce_out_curs() for n in m.snk_nodeblocks()])
    return {'src_in': src_in,
            'src_out': src_out,
            'snk_in': snk_in,
            'snk_out': snk_out}


def aggregate_current_vectors(m):
    ms = m.h
    return (np.array([[m.nodes[j][i].aggregate_current_vector()[0]
                     for i in range(ms)] for j in range(ms)]),
            np.array([[m.nodes[j][i].aggregate_current_vector()[1]
                     for i in range(ms)] for j in range(ms)]))


def nodal_current_dict(node):
    return {d: node.ammeters[d].current if d in node.ammeters else 0.
            for d in node.directions}


def node_currents(m):
    for n in m.iter_nodes():
        print(nodal_current_dict(n))



def plot_virtual_heatmap(m, current_flow_plot=None):
    import matplotlib.pyplot as plt

    ms = m.hp.N
    potentials = m.final_grid
    fig, axes = plt.subplots(1, 1)
    axes.imshow(potentials, cmap='hot', interpolation='nearest')

    plt.show()

def plot_heatmap(m, current_flow_plot=None):
    import matplotlib.pyplot as plt

    ms = m.h
    potentials = m.final_grid 
    fig, axes = plt.subplots(1, 1)
    axes.imshow(potentials, cmap='hot', interpolation='nearest')

    if current_flow_plot is not None:
        U, V = aggregate_current_vectors(m)
        if current_flow_plot == 'stream':
            axes.streamplot(np.array([i for i in range(ms)]),
                            np.array([i for i in range(ms)]),
                            U, -V, color='green', linewidth=1,
                            density=1)
        elif current_flow_plot == 'quiver':
            axes.quiver(np.array([i for i in range(ms)]),
                        np.array([i for i in range(ms)]),
                        U, V, color='green')
        else:
            print('Unrecognized current_flow_plot')

    plt.show()

def plot_surface(m):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm

    potentials = m.final_grid 
    gs = potentials.shape[0]

    fig = plt.figure()
    axes = fig.gca(projection='3d')

    x,y = np.meshgrid([i for i in range(gs)],
                      [i for i in range(gs)])

    axes.plot_surface(x, y, potentials,
                      rstride=1, cstride=1,
                      cmap=cm.coolwarm, linewidth=0,
                      antialiased=True, shade=True)

    plt.show()

def plot_errmap(data1, data2):
    err = abs(data1-data2)

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 1)

    axes.imshow(err, cmap='hot', interpolation='nearest')

    plt.show()

def print_error_table(sim_grid, base):
    gs = sim_grid.shape[0]

    frs = '{0:>9}    {1:>12}    {2:>12}    {3:>12}'
    print(frs.format('Node', 'Potential', 'Base', 'ERROR'))
    for i,j in it.product(range(gs), range(gs)):
        s = sim_grid[i,j]
        b = base[i,j]
        print(frs.format(str((i,j)), '{: e}'.format(s),
                                     '{: e}'.format(b),
                                     '{: e}'.format(s-b)))

# returns a list of lists of dictionaries of dictionaries:
# run_current_split_analysis(model)[i][j]['E'] sis a dictionary that
# stores the current split at node ij when the current is incoming from
# East
# Note: The model doesn't need to have been run prior to calling this
def run_current_split_analysis(m, verbose=False):
    ms = m.h

    m.run_spice_solver()

    # record initial current splits
    init_nodal_current_splits = [[nodal_current_dict(m.nodes[i][j])
                                 for j in range(ms)] for i in range(ms)]

    final_nodal_current_splits = [[{} for j in range(ms)]
                                  for i in range(ms)]

    # print(init_nodal_current_splits)

    def split_dif(biased, base):
        return {d: abs(biased[d]-base[d]) for d in ('E', 'W', 'N', 'S')}

    def mask_split(split, mask_key):
        split[mask_key] = 0.0

    def normalize_split(split):
        val_sum = sum([v for _, v in split.items()])
        if val_sum != 0:
            return {k: v/val_sum for k, v in split.items()}
        else:
            return split

    # apply bias and measure the deltas
    for i, j in it.product(range(ms), range(ms)):
        node = m.nodes[i][j]
        base_split = init_nodal_current_splits[i][j]
        fin_splits = final_nodal_current_splits[i][j]

        for d in node.directions:
            if d in node.ammeters:
                a = node.ammeters[d]
                if verbose:
                    sys.stdout.write(
                            '\rRunning bias: {}, {}'.format(str((i, j)),
                                                            d))
                a.set_bias(1)
                m.run_spice_solver()
                biased_split = nodal_current_dict(node)
                a.reset_bias()
                dif_split = split_dif(biased_split, base_split)
                mask_split(dif_split, d)
                fin_split = normalize_split(dif_split)
                fin_splits[d] = fin_split
            else:
                fin_splits[d] = {'E': 0., 'W': 0., 'N': 0., 'S': 0.}

    return final_nodal_current_splits


# a little bit different then other library functions this take the
# return value from the curren split analysis function
# The rationale is that I want to separate functions that run the model
# and collect datar from the model in the future. In that vein, a
# printer mustn't run any solvers. But at the same time, I don't think
# there is feasible way to store the curretn split data within the
# model. So the data should be stored in the application level even if
# the whole purpose is just to print the whole data
def print_current_splits(splits):
    mesh_size = len(splits)
    row_format = '{}\t\t{}\n\t\t{}\n\t\t{}\n\t\t{}'
    dict_format = 'E: {E:0<.2f}  W: {W:0<.2f}  N: {N:0<.2f}  S: {S:0<.2f}'
    for i, j in it.product(range(mesh_size), range(mesh_size)):
        tmp_splits = splits[i][j]
        print(row_format.format(str((i, j)),
              *[dict_format.format(**d) for _, d in tmp_splits.items()]))

        print()

def load_grid_from_comsol_csv(filename):
    import csv
    with open(filename, newline='\n') as csvfile:
        reading_data = False
        comsol_reader = csv.reader(csvfile, delimiter=',')
        offset = 0.
        count = 0
        for r in comsol_reader:
            if not reading_data:
                if r[0] == '% Nodes':
                    num_nodes = int(r[1])
                    tmp_grid_size = math.sqrt(num_nodes)
                    grid_size = int(tmp_grid_size)
                    if grid_size != tmp_grid_size:
                        print('Grid is not square')
                    grid = np.zeros((grid_size, grid_size))

                elif r[0] == '% X': # this depends on comsol version
                    reading_data = True
            else: # now we are reading the data
                def to_ij(xy):
                    return (int(grid_size*(xy[1]-offset)),
                            int(grid_size*(xy[0]-offset)))

                tmp_x = float(r[0])
                tmp_y = float(r[1])
                tmp_val = float(r[2])
                if offset == 0.:
                    offset = tmp_x

                grid[int(count/grid_size), count%grid_size] = tmp_val
                count += 1

        return grid
