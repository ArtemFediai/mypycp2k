import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

csv_file = '../data/qm9_postprocessing/output_10.csv'
csv_file = '../data/qm9_postprocessing/output_cluster_10.csv'

data_frame = pd.read_csv(csv_file, delimiter=', ')

# todo: add log of the total time. or N**1/3 as x axis.

sns.boxenplot(data=data_frame,
              x='n_electrons',
              y='cpu_sec',
              hue='cpu_model',
              showfliers=False
              )

data_frame[r'log$_{10}$(cpu_sec)'] = np.log10(data_frame['cpu_sec'])
data_frame[r'log$_{10}(cpu hours)$'] = np.log10(data_frame['cpu_sec'] / 360)
data_frame['cpu hours'] = data_frame['cpu_sec'] / 360
data_frame[r'$n_e^3/10^4$'] = data_frame['n_electrons'] ** 3 / 1E4

plt.ylim([0, 60000])
# plt.xticks([10, 20, 60])

plt.show()

print('done')

sns.boxenplot(data=data_frame,
              x='n_electrons',
              y='cpu hours',
              hue='cpu_model',
              showfliers=False
              )
plt.show()

total_cpu_hours = data_frame['cpu hours'].sum()
print(total_cpu_hours)

sns.lineplot(
    data=data_frame,
    x='n_electrons',
    y='cpu hours',
)

plt.show()

# below the Figure will be used for the paper!

plt.figure(figsize=(6,4))

sns.lineplot(
    data=data_frame,
    x=r'$n_e^3/10^4$',
    y='cpu hours',
    markers=True,
    err_style='bars',
    # label='std of actual data',
    linestyle='None'
)
# plt.show()


sns.regplot(
    data=data_frame,
    x=r'$n_e^3/10^4$',
    y='cpu hours',
    scatter=False,
    x_ci='range',
    ci=None,
    line_kws={'color': 'C0', 'linestyle': '--'},
    label='linear regression',
)

# plt.show()


plt.plot(
    [0, 40], [0, 80],
    color='None',
    linestyle='None',
    # label='cpu hours = $n_e^3$'
)

# second ax

plt.legend()

plt.xlim(left=0)
plt.ylim(bottom=0)

ax1 = plt.gca()
ax2 = ax1.twiny()


def cube_root(x):
    return (x * 1E4) ** (1 / 3)


def my_cube(x):
    return x ** 3 / 1E4

xticks_ax1 = my_cube(np.array([10, 40, 50, 60, 70, 74]))
ax1.set_xticks(xticks_ax1)


ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(ax1.get_xticks())
ax2.set_xticklabels([f'{cube_root(x):.0f}' for x in ax1.get_xticks()])


#
# n_el = [10, 20, 80]
# ax2_ticks = n_el  # yes
# ax1.set_xticks(ax2, n_el)

# ax1.set_xticks([f'{my_cube(x):.2f}' for x in ax2_ticks])

# ax2.set_xlim(ax1.get_xlim())
# ax2.set_xticks(ax1.get_xticks())
# ax2.set_xticklabels([f'{cube_root(x):.2f}' for x in ax1.get_xticks()])
# <-- second ax

ax1.set_title(r'$n_e$')

plt.savefig('cpu_hours_vs_n_e', dpi=600,
            bbox_inches='tight'
            )

plt.show()
# this is good!


data_frame['n_electrons'].max()  # 74
data_frame['n_electrons'].min()  # 10