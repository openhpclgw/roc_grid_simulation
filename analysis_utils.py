import numpy as np
import itertools as it
import sys

def print_current_table(m):
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
    mesh_size = m.h
    return (np.array([[m.nodes[j][i].aggregate_current_vector()[0] for i in
        range(mesh_size)] for j in range(mesh_size)]),
            np.array([[m.nodes[j][i].aggregate_current_vector()[1] for i in
        range(mesh_size)] for j in range(mesh_size)]))

def nodal_current_dict(node):
    return {d:node.ammeters[d].current if d in node.ammeters else 0. 
            for d in node.directions}

def node_currents(m):
    for n in m.iter_nodes():
        print(nodal_current_dict(n))

def node_potentials(m):
    mesh_size = m.h
    return  np.array([[m.nodes[j][i].potential for i in
        range(mesh_size)] for j in range(mesh_size)])

def plot_heatmap(m, current_flow_plot=None):
    import matplotlib.pyplot as plt

    mesh_size = m.h
    potentials = node_potentials(m)
    fig, axes = plt.subplots(1,1)
    axes.imshow(potentials, cmap='hot', interpolation='nearest')

    if current_flow_plot is not None:
        U, V = aggregate_current_vectors(m)
        if current_flow_plot == 'stream':
            axes.streamplot(np.array([i for i in range(mesh_size)]),
                           np.array([i for i in range(mesh_size)]),
                           U, -V, color='green', linewidth=1, density=1)
        elif current_flow_plot == 'quiver':
            axes.quiver(np.array([i for i in range(mesh_size)]),
                           np.array([i for i in range(mesh_size)]),
                           U, V, color='green')
        else:
            print('Unrecognized current_flow_plot')

    plt.show()

# returns a list of lists of dictionaries of dictionaries:
# run_current_split_analysis(model)[i][j]['E'] sis a dictionary that
# stores the current split at node ij when the current is incoming from
# East
# Note: The model doesn't need to have been run prior to calling this
def run_current_split_analysis(m):
    mesh_size = m.h

    m.run_spice_solver()

    # record initial current splits
    init_nodal_current_splits = [[nodal_current_dict(m.nodes[i][j])
            for j in range(mesh_size)] for i in range(mesh_size)]

    final_nodal_current_splits = [[{}
            for j in range(mesh_size)] for i in range(mesh_size)]

    # print(init_nodal_current_splits)

    def split_dif(biased, base):
        return {d:abs(biased[d]-base[d]) for d in ('E', 'W', 'N', 'S')}

    def mask_split(split, mask_key):
        split[mask_key] = 0.0

    def normalize_split(split):
        val_sum = sum([v for _,v in split.items()])
        if val_sum != 0:
            return {k:v/val_sum for k,v in split.items()}
        else:
            return split

    # apply bias and measure the deltas
    for i,j in it.product(range(mesh_size), range(mesh_size)):
        node = m.nodes[i][j]
        base_split = init_nodal_current_splits[i][j]
        fin_splits = final_nodal_current_splits[i][j]

        for d in node.directions:
            if d in node.ammeters:
                a = node.ammeters[d]
                sys.stdout.write('\rRunning bias: {}, {}'.format(str((i,j)), d))
                a.set_bias(1)
                m.run_spice_solver()
                biased_split = nodal_current_dict(node)
                a.reset_bias()
                dif_split = split_dif(biased_split, base_split)
                mask_split(dif_split, d)
                fin_split = normalize_split(dif_split)
                fin_splits[d] = fin_split
            else:
                fin_splits[d] = {'E':0.,'W':0.,'N':0.,'S':0.}

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
    dict_format = 'E:{E:0<.2f} W:{W:0<.2f} N:{N:0<.2f} S:{S:0<.2f}'
    for i,j in it.product(range(mesh_size), range(mesh_size)):
        tmp_splits = splits[i][j]
        print(row_format.format(str((i,j)),
            *[dict_format.format(**d) for _,d in tmp_splits.items()]))

        
        print()



