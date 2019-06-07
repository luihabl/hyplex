import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from scipy.constants import m_e, e, epsilon_0

dens_i = np.mean(np.transpose(pd.read_csv('dens_i.csv', header=None).values), axis=0)
dens_e = np.mean(np.transpose(pd.read_csv('dens_e.csv', header=None).values), axis=0)

di = dens_i
de = dens_e
diff = di - de

x = np.arange(dens_i.shape[0])

fig, ax = plt.subplots(nrows=1, figsize=(6, 5))

#ax.plot(x, di, 'b', label='Ions')
#ax.plot(x, de, 'r', label='Electrons')
ax.plot(x, diff, 'b', label='Difference')

ax.xaxis.set_major_locator(plt.MultipleLocator(64))

#plt.axhline(1, color="k", linestyle='--')
plt.ylabel('$\Delta n$ [$m^{-3}$]')
plt.xlabel('x/$\lambda_{De}$')

plt.xlim( (0, 256) )
#plt.tight_layout()
plt.show()
