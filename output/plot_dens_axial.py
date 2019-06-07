import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from scipy.constants import m_e, e, epsilon_0

dens_i = pd.read_csv('dens_i.csv', header=None).values
dens_e = pd.read_csv('dens_e.csv', header=None).values

di = dens_i[:, 0]
de = dens_e[:, 0]

x = np.arange(dens_i.shape[0])

fig, ax = plt.subplots(nrows=1, figsize=(6, 5))

ax.semilogy(x, di, 'b', label='Ions')
ax.semilogy(x, de, 'r', label='Electron')

ax.xaxis.set_major_locator(plt.MultipleLocator(64))

plt.axhline(1e16, color="k", linestyle='--')
plt.ylabel('$n$ [$10^{16}$ $m^{-3}$]')
plt.xlabel('x/$\lambda_{De}$')

plt.xlim( (0, 256) )
plt.ylim( (1e14, 1e17) )
#plt.tight_layout()
plt.legend()
plt.show()
