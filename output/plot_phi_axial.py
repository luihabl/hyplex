import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from scipy.constants import m_e, e, epsilon_0

phi = np.transpose(pd.read_csv('phi.csv', header=None).values)

p = phi[0, :]

x = np.arange(phi.shape[1])

fig, ax = plt.subplots(nrows=1, figsize=(6, 5))

ax.plot(x, p, 'b')

ax.xaxis.set_major_locator(plt.MultipleLocator(64))
ax.yaxis.set_major_locator(plt.MultipleLocator(10))

plt.axhline(0, color="k", linestyle='--')
plt.ylabel('$\phi$ [V]')
plt.xlabel('x/$\lambda_{De}$')
plt.ylim( (-20, 20) )
plt.xlim( (0, 256) )
#plt.tight_layout()
plt.show()
