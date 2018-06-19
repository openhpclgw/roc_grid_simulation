import matplotlib.pyplot as plt
import math
import itertools as it
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--log-scale-x', action='store_true')
args = parser.parse_args()

x = (10, 10000, 100000, 1000000)

cases=('1', '2', '3')
sizes=('8', '16', '32')
err_file='case{case}/size{size}/errors'

title='Case:{case}, Size: {size}x{size}'

fig, axs = plt.subplots(ncols=len(sizes), nrows=len(cases),
                        sharex=True, sharey=True,
                        figsize=(12,8))

plt.subplots_adjust(top=0.85, bottom=0.15, hspace=0.3)

base_ax = fig.add_subplot(111, frameon=False)

get_lines = True
lines = []
for i, (c,s) in enumerate(it.product(cases, sizes)):
    print(i)
    ax_i, ax_j = i//len(sizes), i%len(sizes)
    ax = axs[ax_i, ax_j]

    maxs = []
    means = []
    medians = []

    with open(err_file.format(case=c, size=s)) as f:
        for line in f:
            for val, lst in zip(line.split(), (maxs, means, medians)):
                lst.append(float(val))

    if args.log_scale_x:
        ax.set_xscale('log')
    ax.set_title(title.format(case=c, size=s), fontsize=16)
    line1, = ax.plot(x, maxs, label='Max')
    line2, = ax.plot(x, means, label='Mean')
    line3, = ax.plot(x, medians, label='Median')
    if get_lines:
        lines.append(line1)
        lines.append(line2)
        lines.append(line3)
        get_lines = False
    ax.set_xticks(x)
    if args.log_scale_x:
        ax.set_xticklabels((r'$10^0$', r'$10^3$', r'$10^4$', r'$10^5$'))
    else:
        ax.set_xticklabels(('', '', r'$10^4$', r'$10^5$'))
    ax.set_ylim((0.,1.))

# plt.tick_params(labelcolor='none', top='off', bottom='off', left='off',
        # right='off')
# plt.grid(False)
base_ax.set_xticks(())
base_ax.set_yticks(())
base_ax.set_xticklabels(labels=[''])
base_ax.set_xlabel(r'Loss Multiplier (dB/$\Omega$)', labelpad=50)
base_ax.set_ylabel('Error', labelpad=50)


plt.figlegend(handles=lines, labels=('Max', 'Mean', 'Median'),
              ncol=3, loc='upper center', fontsize=16)
# plt.tight_layout()
plt.show()

