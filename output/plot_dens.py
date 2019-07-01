import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from scipy.constants import m_e, e, epsilon_0
import scipy.ndimage
m_i = 6.7e-27


dens_i = np.transpose(pd.read_csv('dens_i.csv', header=None).values)
dens_e = np.transpose(pd.read_csv('dens_e.csv', header=None).values)

x = np.arange(dens_e.shape[1])
y = np.arange(dens_e.shape[0])
x, y = np.meshgrid(x, y)

fig, axes = plt.subplots(nrows=2, figsize=(10, 6))

i = axes[0].imshow(dens_i, cmap='plasma', origin='lower')
e = axes[1].imshow(dens_e, cmap='plasma', origin='lower')

axes[0].yaxis.set_major_locator(plt.MultipleLocator(16))
axes[1].yaxis.set_major_locator(plt.MultipleLocator(16))

axes[0].xaxis.set_major_locator(plt.MultipleLocator(64))
axes[1].xaxis.set_major_locator(plt.MultipleLocator(64))

divider0 = make_axes_locatable(axes[0])
divider1 = make_axes_locatable(axes[1])
cax0 = divider0.append_axes("right", size="1%", pad=0.1)
cax1 = divider1.append_axes("right", size="1%", pad=0.1)

axes[0].axis('image')
axes[0].set(title='Ion density')
axes[0].set_ylabel('y/$\lambda_{De}$')
fig.colorbar(i, ax=axes[0],cax=cax0, label='$n_i$ [$m^{-3}$]')

axes[1].axis('image')
axes[1].set(title='Electron density')
axes[1].set_ylabel('y/$\lambda_{De}$')
axes[1].set_xlabel('x/$\lambda_{De}$')
fig.colorbar(e, ax=axes[1],cax=cax1, label='$n_e$ [$m^{-3}$]')

plt.tight_layout()
plt.show()
