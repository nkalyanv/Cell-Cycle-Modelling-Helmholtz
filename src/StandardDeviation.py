from tracemalloc import start
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

YPD = np.load('YPDprotein.npy')
SCD = np.load('SCDprotein.npy')
SCGE = np.load('SCGEprotein.npy')

YPD_mean = np.mean(YPD)
SCD_mean = np.mean(SCD)
SCGE_mean = np.mean(SCGE)

YPD_std = np.std(YPD)
SCD_std = np.std(SCD)
SCGE_std = np.std(SCGE)

mediums = ['SCGE', 'SCD', 'YPD']
x_pos = np.arange(len(mediums))
CTEs = [SCGE_mean, SCD_mean, YPD_mean]
error = [SCGE_std, SCD_std, YPD_std]


fig, ax = plt.subplots()
ax.bar(x_pos, CTEs, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel('Protein Content')
ax.set_xticks(x_pos)
ax.set_xticklabels(mediums)
#ax.yaxis.grid(True)

# Save the figure and show
plt.tight_layout()
plt.savefig('Plots/proteinContent_barPlot.png')