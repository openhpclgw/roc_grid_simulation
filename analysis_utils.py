import numpy as np

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

