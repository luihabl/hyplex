import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib

import h5py
import numpy as np
import os

from scipy.constants import m_e, e, epsilon_0, k
m_i = 2.17e-25
at_sccm = 4.477962e17

def plot_field(field, label=''):
    fig, axes = plt.subplots(nrows=1, figsize=(10, 4))

    p = axes.imshow(field, cmap='plasma', origin='lower')

    axes.xaxis.set_major_locator(plt.MultipleLocator(64))
    axes.yaxis.set_major_locator(plt.MultipleLocator(16))

    divider0 = make_axes_locatable(axes)
    cax0 = divider0.append_axes("right", size="1%", pad=0.1)
    
    fig.colorbar(p, ax=axes,cax=cax0, label=label)
    
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


def load_dataset(h5_file, key):
    try:
        data = h5_file[key][()]
        if data.shape[1] > 1:
            data = data.transpose()
        return data
    except:
        print('failed to load ' + h5_file.filename + ' / ' + key)
        return None


def load_file_data(file_path):
    d = {}
    file = h5py.File(file_path)
    for k in file.keys():
        d[k] = load_dataset(file, k)
    return d
    
def load_all(folder_path='.'):
    d = {}
    for f in os.listdir(folder_path):
        if f.endswith('.h5'):
            print('Loaded ' + f)
            d[f.split('.')[0]] = load_file_data(os.path.join(folder_path, f))
    return d


if __name__=='__main__':
    global d
    d = load_all()
