import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib

import pandas as pd
import numpy as np

from scipy.constants import m_e, e, epsilon_0, k
m_i = 2.17e-25
at_sccm = 4.477962e17

def plot_field(field):
    fig, axes = plt.subplots(nrows=1, figsize=(10, 4))

    p = axes.imshow(field, cmap='plasma', origin='lower')

    axes.xaxis.set_major_locator(plt.MultipleLocator(64))
    axes.yaxis.set_major_locator(plt.MultipleLocator(16))

    divider0 = make_axes_locatable(axes)
    cax0 = divider0.append_axes("right", size="1%", pad=0.1)
    
    fig.colorbar(p, ax=axes,cax=cax0)
    
    axes.axis('image')
    axes.set_ylabel('y/$\lambda_{De}$')
    axes.set_xlabel('x/$\lambda_{De}$')

    plt.tight_layout()
    plt.show()
    
def plot_scatter(x, y, color='white', background='dark'):
    if(background == 'dark'):
        with plt.style.context('dark_background'):
            fig, axes = plt.subplots(nrows=1, figsize=(12, 3)) 
            axes.plot(x, y,color=color,marker=',',lw=0, linestyle="")
    else:
        fig, axes = plt.subplots(nrows=1, figsize=(12, 3)) 
        axes.plot(x, y,color=color,marker=',',lw=0, linestyle="")
    plt.tight_layout()
    plt.show()
    
def plot_distribution(v, bins='auto'):
    plt.hist(v, bins=bins, histtype='step', color='r')
    plt.show()

# def plot_line(x=None, y):
#     x = np.arange(phi.shape[0]) if x == None else x
#     plt.plot(x, y)
#     plt.show()


def load_fmatrix(path):
    try:
        return np.transpose(pd.read_csv(path, header=None).values)
    except:
        print('failed to load ' + path)
        return None

def load_array(path):
    try:
        return pd.read_csv(path, header=None).values
    except:
        print('failed to load ' + path)
        return None

def load():
    global dens_e, dens_i, dens_n, p_e, p_i, phi, phi_av, v_e, v_i, k_e, k_i, ne, ni, v_cap
    
    dens_e = load_fmatrix('dens_e_state.csv')
    dens_i = load_fmatrix('dens_i_state.csv')
    dens_n = load_fmatrix('dens_n.csv')
    phi = load_fmatrix('phi_state.csv')
    phi_av = load_fmatrix('phi_av.csv')
    ne = load_array('n_active_e.csv')
    ni = load_array('n_active_i.csv')
    v_cap = load_array('v_cap.csv')
    p_e = load_fmatrix('p_e_state.csv')
    p_i = load_fmatrix('p_i_state.csv')

    try:    
        v_e = p_e[3:6, :]
        k_e = 0.5 * m_e * (v_e[0] ** 2 + v_e[1] ** 2 + v_e[2] ** 2)
    except:
        print('failed to calc v_e and k_e')

    try:
        v_i = p_i[3:6, :]
        k_i = 0.5 * m_i * (v_i[0] ** 2 + v_i[1] ** 2 + v_i[2] ** 2)
    except:
        print('failed to calc v_i and k_i')

    

if __name__=='__main__':
    load()
    x = np.arange(phi.shape[1]) * 0.166e-3
    y = np.arange(phi.shape[0]) * 0.166e-3
    X, Y = np.meshgrid(x, y)

