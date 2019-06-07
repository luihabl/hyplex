import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from scipy.constants import m_e, e, epsilon_0

dens_i = pd.read_csv('dens_i.csv', header=None).values



di = dens_i[15, :]
y = np.arange(2 * dens_i.shape[1]) - dens_i.shape[1]

fig, ax = plt.subplots(nrows=1, figsize=(6, 5))


K = 5
for i in range(K):
    ax.plot(y, np.append(np.flip(dens_i[i * int(256 / (K- 1)),:]), dens_i[i * int(256 / (K- 1)), :]), label= str(i * int(256 / (K- 1))))

ax.xaxis.set_major_locator(plt.MultipleLocator(64))

plt.ylabel('$n$ [ $m^{-3}$]')
plt.xlabel('y/$\lambda_{De}$')

#plt.tight_layout()
plt.legend()
plt.show()
