import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from scipy.constants import m_e, e, epsilon_0

phi = np.transpose(pd.read_csv('phi.csv', header=None).values)

p = phi[0, :]

x = np.arange(phi.shape[1]) * 0.166e-3

phi_t = (-2e14 * e / epsilon_0) * (0.5 * x ** 2 - 42.49e-3 * x) 

fig, ax = plt.subplots(nrows=1, figsize=(6, 5))

ax.plot(x, p, 'b')
ax.plot(x, phi_t, 'r')

ax.xaxis.set_major_locator(plt.MultipleLocator(64))

plt.axhline(0, color="k", linestyle='--')
plt.ylabel('$\phi$ [V]')
plt.xlabel('x/$\lambda_{De}$')
#plt.ylim( (-20, 20) )
#plt.xlim( (0, 256) )
#plt.tight_layout()
plt.show()
