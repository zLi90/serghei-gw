"""
    Plot results for the discharge at the outlet
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# output file directory
fdir = ['output-norain-test-2', 'output-rain']
lgd = ['ParFlow', 'CATHY', 'HGS', 'Cast3M', 'SERGHEI']
dataname = 'discharge/'

# Read data
data_files = {
    'pf': ['PF-norain.csv', 'PF-rain.csv'],
    'cy': ['CATHY-norain.csv', 'CATHY-rain.csv'],
    'hgs': ['HGS-norain.csv', 'HGS-rain.csv'],
    'cast1': ['cast3m-norain.csv', None]
}

data_dict = {key: [pd.read_csv(dataname + fname, header=None) if fname else None for fname in fnames]
             for key, fnames in data_files.items()}

# Read domainTimeSeries data
data = [np.genfromtxt(fdir[ff] + '/domainTimeSeries.out', delimiter=' ', skip_header=True) for ff in range(len(fdir))]

# Colors and markers
colors = [(19/255, 103/255, 158/255), (171/255, 58/255, 41/255), (208/255, 127/255, 44/255), (111/255, 109/255, 161/255)]
markers = ['o', '^', 's', '*']
titles = ["Scenario 1", "Scenario 2"]

# Plotting
plt.figure(1, figsize=[12, 4])

for i in range(2):
    plt.subplot(1, 2, i+1)
    for j, (key, marker) in enumerate(zip(data_files.keys(), markers)):
        if data_dict[key][i] is not None:
            plt.scatter(data_dict[key][i].iloc[:, 0], data_dict[key][i].iloc[:, 1], s=14, marker=marker, facecolor='None', edgecolor=colors[j])
    
    plt.plot(data[i][:, 0]/3600, 2.0 * data[i][:, 4] * 3600, color='k')
    plt.title(titles[i])
    plt.xlabel('Time [h]')
    plt.ylabel('Flow Rate [$m^{3}/h$]')
    plt.xlim([0, 120])
    plt.ylim([0, 10] if i == 0 else [0, 1200])

    if i==0:
        plt.legend(lgd)

plt.savefig('discharge.png', format='png', bbox_inches='tight', dpi=600)
plt.show()
