import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from scipy.constants import m_e, e, epsilon_0

dens_i = np.transpose(pd.read_csv('dens_i.csv', header=None).values)
dens_e = np.mean(np.transpose(pd.read_csv('dens_e.csv', header=None).values), axis=1)

d = dens_i[0, :]

x = np.arange(dens_i.shape[1])

fig, ax = plt.subplots(nrows=1, figsize=(6, 5))

ax.semilogy(x, d/1e16, 'b')

ax.xaxis.set_major_locator(plt.MultipleLocator(64))

plt.axhline(1, color="k", linestyle='--')
plt.ylabel('$n_i$ [$10^{16}$ $m^{-3}$]')
plt.xlabel('x/$\lambda_{De}$')
plt.ylim( (1e-2, 1e1) )
plt.xlim( (0, 256) )
#plt.tight_layout()
plt.show()
