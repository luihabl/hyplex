import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib

import h5py
import exdir
import numpy as np
import argparse

from functools import reduce
import operator
from pathlib import Path

from scipy.constants import m_e, e, epsilon_0, k

m_i = 2.17e-25
at_sccm = 4.477962e17


class OutputFile:
    
    def __init__(self, _file, _sep='.'):
        self.file = _file
        self.sep = _sep

    def __getitem__(self, path):
        return reduce(operator.getitem, path.split(self.sep), self.file)

    def field_img(self, path, label=''):

        field = np.transpose(self[path])

        fig, axes = plt.subplots(nrows=1, figsize=(10, 4))

        p = axes.imshow(field, cmap='plasma', origin='lower')

        axes.xaxis.set_major_locator(plt.MultipleLocator(64))
        axes.yaxis.set_major_locator(plt.MultipleLocator(16))

        divider0 = make_axes_locatable(axes)
        cax0 = divider0.append_axes("right", size="1%", pad=0.1)
        
        fig.colorbar(p, ax=axes,cax=cax0, label=label)
        
        axes.axis('image')
        axes.set_ylabel('$N_y$')
        axes.set_xlabel('$N_x$')

        plt.tight_layout()
        plt.show()

    def scatter(self, path, cols=(0,1), setlim=True, color='white', background='dark'):
        x = self[path][:,cols[0]]
        y = self[path][:,cols[1]]

        if(background == 'dark'):
            with plt.style.context('dark_background'):
                fig, axes = plt.subplots(nrows=1, figsize=(12, 3)) 
                axes.plot(x, y,color=color,marker=',',lw=0, linestyle="")
        else:
            fig, axes = plt.subplots(nrows=1, figsize=(12, 3)) 
            axes.plot(x, y,color=color,marker=',',lw=0, linestyle="")
        
        if setlim:
            plt.xlim((self['mesh.x'].value.min(), self['mesh.x'].value.max()))
            plt.ylim((self['mesh.y'].value.min(), self['mesh.y'].value.max()))
        plt.tight_layout()
        plt.show()
    
    def field_contour(self, path, num=50, label=''):
        field = np.transpose(self[path])
        x = np.transpose(self['mesh.x'])
        y = np.transpose(self['mesh.y'])
        
        fig, axes = plt.subplots(nrows=1, figsize=(10, 4))
        p = axes.contourf(x,y, field, num, cmap='plasma')

        divider0 = make_axes_locatable(axes)
        cax0 = divider0.append_axes("right", size="1%", pad=0.1)

        fig.colorbar(p, ax=axes,cax=cax0, label=label)

        axes.axis('image')
        axes.set_ylabel('$y (m)$')
        axes.set_xlabel('$x (m)$')

        plt.show()

    def plot_line(self, x1_path, x2_path=None):
        if(x2_path is None):
            plt.plot(self[x1_path].value, 'r')
        else:
            plt.plot(self[x1_path].value, self[x2_path].value, 'r')
        plt.show()
    
    def dist(self, path, col=0, bins='auto'):
        plt.hist(self[path][:,col], bins=bins, histtype='step', color='r')
        plt.show()

    def avg_ke(self, path, mass, cols=[3, 4, 5]):
        t = np.zeros(len(cols))
        for i, c in enumerate(cols):
            t[i] = np.mean(0.5 * self[path][:,c]**2 * mass) / e
        return t


def field_img(field, label=''):

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

def scatter(x, y, color='white', background='dark'):
    if(background == 'dark'):
        with plt.style.context('dark_background'):
            fig, axes = plt.subplots(nrows=1, figsize=(12, 3)) 
            axes.plot(x, y,color=color,marker=',',lw=0, linestyle="")
    else:
        fig, axes = plt.subplots(nrows=1, figsize=(12, 3)) 
        axes.plot(x, y,color=color,marker=',',lw=0, linestyle="")
    plt.tight_layout()
    plt.show()
    
def dist(v, bins='auto'):
    plt.hist(v, bins=bins, histtype='step', color='r')
    plt.show()

def load_file_data_h5(file_path):
    return h5py.File(file_path)
    
def load_file_data_exdir(file_path):
    return exdir.File(file_path)

def load_all(folder_path='.'):
    d = {}
    folder = Path(folder_path)
    for f in folder.glob('*'):
        if f.suffix == '.h5':
            print(f'Loaded {f}')
            d[f.stem] = OutputFile(load_file_data_h5(f))
        if f.suffix == '.exdir':
            print(f'Loaded {f}')
            d[f.stem] = OutputFile(load_file_data_exdir(f))
    
    return d


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Loads all the exdir or hdf5 files in a dir and provides plotting functions.')
    parser.add_argument("-f", "--folder", help="directory with the exdir or hdf5 files.", default=".")
    args = parser.parse_args()

    print(args.folder)

    global d
    d = load_all(args.folder)
