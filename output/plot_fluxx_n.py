import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from scipy.constants import m_e, e, epsilon_0, k

flux_x = np.transpose(pd.read_csv('fluxn_x.csv', header=None).values)
dens = np.transpose(pd.read_csv('dens_n.csv', header=None).values)
wmesh = np.transpose(pd.read_csv('wmesh_n.csv', header=None).values)
vmesh = np.transpose(pd.read_csv('vmesh.csv', header=None).values)
#dens = dens * k * 300 * 7.5
x = np.arange(flux_x.shape[1])
y = np.arange(flux_x.shape[0])
x, y = np.meshgrid(x, y)

fig, axes = plt.subplots(nrows=1, figsize=(10, 4))

p = axes.imshow(flux_x/dens, cmap='plasma', origin='lower')

axes.xaxis.set_major_locator(plt.MultipleLocator(64))
axes.yaxis.set_major_locator(plt.MultipleLocator(16))

divider0 = make_axes_locatable(axes)
cax0 = divider0.append_axes("right", size="1%", pad=0.1)

axes.axis('image')
axes.set(title='Velocity field in x direction')
axes.set_ylabel('y/$\lambda_{De}$')
axes.set_xlabel('x/$\lambda_{De}$')
#fig.colorbar(p, ax=axes,cax=cax0, label='$P$ [mTorr]')
fig.colorbar(p, ax=axes,cax=cax0, label='$u_x$ [$m s^{-1}$]')

plt.tight_layout()
plt.show()
