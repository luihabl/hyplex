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

    def field_img(self, path, label='', inner_label='', save=False, name='field.pdf', multiplier=1, add_factor=0, **kwargs):

        field = np.transpose(self[path])

        fig, axes = plt.subplots(nrows=1, figsize=(10, 4))

        p = axes.imshow(field * multiplier + add_factor, interpolation='nearest', origin='lower', **kwargs)
        
#         axes.xaxis.set_major_locator(plt.MultipleLocator(0.5))
#         axes.yaxis.set_major_locator(plt.MultipleLocator(0.5))
#         axes.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
#         axes.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

#         divider0 = make_axes_locatable(axes)
#         cax0 = divider0.append_axes("right", size="1%", pad=0.1)
#         fig.colorbar(p, ax=axes,cax=cax0, label=label)
        fig.colorbar(p, fraction=0.01225, pad=0.01, label=label)
    
#         axes.axis('image')
        axes.text(0.98, 0.78, inner_label, horizontalalignment='right', verticalalignment='bottom', transform=axes.transAxes, color='white')
        axes.set_ylabel('$y / L_y$')
        axes.set_xlabel('$x / L_x$')
        
        if save:
            fig.savefig(name,  bbox_inches='tight')

        plt.tight_layout()
        plt.show()
        
        return fig, axes

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
    
    def field_contour(self, path, num=50, label='', **kwargs):
        field = np.transpose(self[path])
#         x = np.transpose(self['mesh.x'])
#         y = np.transpose(self['mesh.y'])
        x = np.arange(512)
        y = np.arange(128)
        
        fig, axes = plt.subplots(nrows=1, figsize=(10, 4))
        p = axes.contourf(x,y, field, num, **kwargs)

        divider0 = make_axes_locatable(axes)
        cax0 = divider0.append_axes("right", size="1%", pad=0.1)

        fig.colorbar(p, ax=axes,cax=cax0, label=label)

        axes.axis('image')
        axes.set_ylabel('$y / \lambda_D$')
        axes.set_xlabel('$x / \lambda_D$')
        plt.show()

    def plot_line(self, x1_path, x2_path=None):
        if(x2_path is None):
            plt.plot(self[x1_path].value, 'r')
        else:
            plt.plot(self[x1_path].value, self[x2_path].value, 'r')
        plt.show()
    
    def plot_series(self, y_path, xlabel=None, ylabel=None, xmag=None, title=None, **kwargs):
        path_comp = y_path.split(self.sep)
        t = self[self.sep.join(path_comp[:-1])]['time'].value
        y = self[y_path].value
        
        if xlabel is None:
            plt.xlabel('Time [s]')
        else:
            plt.xlabel(xlabel)
        
        if title is not None:
            plt.title(title)
        
        if ylabel is not None:
            plt.ylabel(ylabel)
        
        tdiv = 1
        if xmag is not None:
            tdiv = xmag
        
        plt.plot(t / tdiv, y, **kwargs)
        plt.tight_layout()
        
        
    
    def dist(self, path, col=0, bins='auto'):
        plt.hist(self[path][:,col], bins=bins, histtype='step', color='r')
        plt.show()

    def avg_ke(self, path, mass, cols=[3, 4, 5]):
        t = np.zeros(len(cols))
        for i, c in enumerate(cols):
            t[i] = np.mean(0.5 * self[path][:,c]**2 * mass) / e
        return t
        
    def plot_vdist(self, vdist_path, title=None, **kwargs):
        path_comp = vdist_path.split(self.sep)
        vlim = self[self.sep.join(path_comp[:-1])].attrs['vlim_' + path_comp[-1]]
        n_v = self[self.sep.join(path_comp[:-1])].attrs['n_v']
        v = np.linspace(vlim[0], vlim[1], n_v)
        
        vdist = self[vdist_path].value
        
        plt.ylabel('Norm. distribution')
        plt.xlabel('v [m/s]')
        
        if title is not None:
            plt.title(title)
        
        plt.plot(v, vdist / vdist.max(), **kwargs)
        
        plt.legend()
        plt.tight_layout()
        
        
def gather_data(d, key):
    v_list = []
    for k in sorted(d):
        v = d[k][key][()].copy()
        v_list.append(v)
    return np.concatenate(v_list)

def n(a, n_avg):
    return np.mean(a.reshape(-1, n_avg), axis=1)


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


def load_latest(folder_path='.'):
    df = None
    folder = Path(folder_path)
    file_list = sorted(folder.glob('*'))
    f = file_list[-1]
    if f.suffix == '.exdir':
        print(file_list[-1])
        return OutputFile(load_file_data_exdir(file_list[-1]))
    else:
        print(f'Incompatible format: {f}')


def load_list(base_path, name_list):
    base_folder = Path(base_path)    
    d = []
    for n in name_list:
        d.append(load_latest(base_folder / n))
    return d
        
if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Loads all the exdir or hdf5 files in a dir and provides plotting functions.')
    parser.add_argument("-f", "--folder", help="directory with the exdir or hdf5 files.", default=".")
    args = parser.parse_args()

    print(args.folder)

    global d
    d = load_all(args.folder)
