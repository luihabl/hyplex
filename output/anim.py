import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

import pandas as pd
import numpy as np
import h5py
import os, sys

from scipy.constants import m_e, e, epsilon_0, k

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
    
        
    
def load_dataset(h5_file, key):
    try:
        data = h5_file[key][()]
        if data.shape[1] > 1:
            data = data.transpose()
        return data
    except:
        print('failed to load ' + h5_file.filename + ' / ' + key)
        return None
        
def load_attribute(h5_file, key, attr_name):
    return h5_file[key].attrs[attr_name]


def load_file_data(file_path):
    d = []
    t = []
    file = h5py.File(file_path, 'r')
    for k in file.keys():
        d.append(load_dataset(file, k))
        t.append(load_attribute(file, k, 'Time [s]')[0])
    return d, t
    
def load_file_data2(file_path):
    d = {}
    file = h5py.File(file_path, 'r')
    for k in file.keys():
        d[k] = load_dataset(file, k)
    return d


def animate(frame, image, axes, data, t):
    print(frame)
    image.set_array(data[frame])
    axes.set_title('Time: {0:.2e} s'.format(t[frame] - t[0]))
    return image

def animate_line(frame, line, axes, data, t):
    print(frame)
    line.set_ydata(data[frame][0:5,:].mean(axis=0))
    axes.set_title('Time: {0:.2e} s'.format(t[frame] - t[0]))
    return (line,)
       
    
def make_line_anim(file_path):
    data, t = load_file_data(file_path)
    subsample = 10
    f0 = 0
    ff = len(data)
    data = data[0:ff//2:subsample]
    t = t[0:ff//2:subsample]
    
    fig, axes = plt.subplots(nrows=1)
    line, = axes.plot(data[0][0:5,:].mean(axis=0))
    # image = axes.imshow(all_data[0], vmin=0, vmax=1)
    vmin = 1e100
    vmax = 0
    
    for d in data:
        if vmax < d[0:5,:].mean(axis=0).max().max(): vmax = d[0:5,:].mean(axis=0).max()
        if vmin > d[0:5,:].mean(axis=0).min().min(): vmin = d[0:5,:].mean(axis=0).min()

    axes.set_ylim(vmin, vmax)
    axes.set_title('Time: {0:.2e} s'.format(t[0] - t[0]))
    
    animation = FuncAnimation(
        # Your Matplotlib Figure object
        fig,
        # The function that does the updating of the Figure
        animate_line,
        # Frame information (here just frame number)
        np.arange(len(data)),
        # Extra arguments to the animate function
        fargs=[line, axes, data, t],
        # Frame-time in ms; i.e. for a given frame-rate x, 1000/x
        interval=1000 / 60
    )
    
    animation.save(file_path.split('.')[0] + '_line' + '.mp4', dpi=150)
    
    
# def make_line_anim(data, t, filename):
#     subsample = 50
#     f0 = 0
#     ff = len(data)
#     data = data[0:ff//2:subsample]
#     t = t[0:ff//2:subsample]
#     
#     fig, axes = plt.subplots(nrows=1)
#     line, = axes.plot(data[0][0:5,:].mean(axis=0))
#     # image = axes.imshow(all_data[0], vmin=0, vmax=1)
#     vmin = 1e100
#     vmax = 0
#     
#     for d in data:
#         if vmax < d[0:5,:].mean(axis=0).max().max(): vmax = d[0:5,:].mean(axis=0).max()
#         if vmin > d[0:5,:].mean(axis=0).min().min(): vmin = d[0:5,:].mean(axis=0).min()
# 
#     axes.set_ylim(vmin, vmax)
#     axes.set_title('Time: {0:.2e} s'.format(t[0] - t[0]))
#     
#     animation = FuncAnimation(
#         # Your Matplotlib Figure object
#         fig,
#         # The function that does the updating of the Figure
#         animate_line,
#         # Frame information (here just frame number)
#         np.arange(len(data)),
#         # Extra arguments to the animate function
#         fargs=[line, axes, data, t],
#         # Frame-time in ms; i.e. for a given frame-rate x, 1000/x
#         interval=1000 / 60
#     )
#     
#     animation.save(filename + '.mp4', dpi=150)

    
def make_anim(file_path):
    
    data, t = load_file_data(file_path)
    subsample = 50
    f0 = 0
    ff = len(data)
    data = data[0:ff//2:subsample]
    t = t[0:ff//2:subsample]
    
    fig, axes = plt.subplots(nrows=1, figsize=(10, 4))
    # image = axes.imshow(all_data[0], vmin=0, vmax=1)
    vmin = 1e100
    vmax = 0
    
    for d in data:
        if vmax < d.max().max(): vmax = d.max()
        if vmin > d.min().min(): vmin = d.min() 
        
    image = axes.imshow(data[0], cmap='plasma', origin='lower', vmin=vmin, vmax=vmax, norm=colors.PowerNorm(gamma=0.5))
    axes.set_title('Time: {0:.2e} s'.format(t[0] - t[0]))
    
    divider0 = make_axes_locatable(axes)
    cax0 = divider0.append_axes("right", size="1%", pad=0.1)
    fig.colorbar(image, ax=axes, cax=cax0)
    
    
    animation = FuncAnimation(
        # Your Matplotlib Figure object
        fig,
        # The function that does the updating of the Figure
        animate,
        # Frame information (here just frame number)
        np.arange(len(data)),
        # Extra arguments to the animate function
        fargs=[image, axes, data, t],
        # Frame-time in ms; i.e. for a given frame-rate x, 1000/x
        interval=1000 / 30
    )
    
    animation.save(file_path.split('.')[0] + '.mp4', dpi=300)


    
# def make_anim(data, t, filename):
#     
#     data = data[::10]
#     t = t[::10]
#     
#     fig, axes = plt.subplots(nrows=1, figsize=(10, 4))
#     # image = axes.imshow(all_data[0], vmin=0, vmax=1)
#     vmin = 1e100
#     vmax = 0
#     
#     for d in data:
#         if vmax < d.max().max(): vmax = d.max()
#         if vmin > d.min().min(): vmin = d.min() 
#         
#     
#     image = axes.imshow(data[0], cmap='plasma', origin='lower', vmin=vmin, vmax=vmax)
#     axes.set_title('Time: {0:.2e} s'.format(t[0] - t[0]))
#     
#     divider0 = make_axes_locatable(axes)
#     cax0 = divider0.append_axes("right", size="1%", pad=0.1)
#     fig.colorbar(image, ax=axes, cax=cax0)
#     
#     
#     animation = FuncAnimation(
#         # Your Matplotlib Figure object
#         fig,
#         # The function that does the updating of the Figure
#         animate,
#         # Frame information (here just frame number)
#         np.arange(1000),
#         # Extra arguments to the animate function
#         fargs=[image, axes, data, t],
#         # Frame-time in ms; i.e. for a given frame-rate x, 1000/x
#         interval=1000 / 30
#     )
#     
#     animation.save(filename + '.mp4', dpi=300)
    


