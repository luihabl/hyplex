import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from scipy.constants import m_e, e, epsilon_0

phi = np.transpose(pd.read_csv('phi.csv', header=None).values)

x = np.arange(phi.shape[1])
y = np.arange(phi.shape[0])
x, y = np.meshgrid(x, y)

fig, axes = plt.subplots(nrows=1, figsize=(10, 4))

p = axes.contourf(x, y, phi, levels=50)

axes.xaxis.set_major_locator(plt.MultipleLocator(64))
axes.yaxis.set_major_locator(plt.MultipleLocator(16))

divider0 = make_axes_locatable(axes)
cax0 = divider0.append_axes("right", size="1%", pad=0.1)

axes.axis('image')
axes.set(title='Potential')
axes.set_ylabel('y/$\lambda_{De}$')
axes.set_xlabel('x/$\lambda_{De}$')
fig.colorbar(p, ax=axes,cax=cax0, label='$\phi$ [V]')


plt.tight_layout()
plt.show()
