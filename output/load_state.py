import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib

import pandas as pd
import numpy as np
import os

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

def load(path='.'):

    data_dict = {}
    
    data_dict['dens_e'] = load_fmatrix(os.path.join(path, 'dens_e_state.csv'))
    data_dict['dens_i'] = load_fmatrix(os.path.join(path, 'dens_i_state.csv'))
    data_dict['dens_n'] = load_fmatrix(os.path.join(path, 'dens_n.csv'))
    data_dict['phi_in'] = load_fmatrix(os.path.join(path, 'phi_state.csv'))
    data_dict['phi_av'] = load_fmatrix(os.path.join(path, 'phi_av.csv'))
    data_dict['ne'] = load_array(os.path.join(path, 'n_active_e.csv'))
    data_dict['ni'] = load_array(os.path.join(path, 'n_active_i.csv'))
    data_dict['v_cap'] = load_array(os.path.join(path, 'v_cap.csv'))
    data_dict['p_e'] = load_fmatrix(os.path.join(path, 'p_e_state.csv'))
    data_dict['p_i'] = load_fmatrix(os.path.join(path, 'p_i_state.csv'))
    
    try:    
        data_dict['k_e'] = 0.5 * m_e * (data_dict['p_e'][3] ** 2 + data_dict['p_e'][4] ** 2 + data_dict['p_e'][5] ** 2) / e
    except:
        print('failed to calc k_e')

    try:
        data_dict['k_i'] = 0.5 * m_i * (data_dict['p_i'][3] ** 2 + data_dict['p_i'][4] ** 2 + data_dict['p_i'][5] ** 2) / e
    except:
        print('failed to calc k_i')
        
    try:
        data_dict['v_cap_av'] = data_dict['v_cap'][-int(0.1 * data_dict['v_cap'].shape[0]):].mean()
    except:
        print('failed to calc k_i')
    
    return data_dict
    

if __name__=='__main__':
    global d
    d = load()
