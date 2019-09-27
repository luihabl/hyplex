import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

import pandas as pd
import numpy as np
import os

from scipy.constants import m_e, e, epsilon_0, k

anim_folder = r'anim/'


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
    
def animate(frame):
    global data, image
    image.set_array(data[frame])
    return image

    
def load_files():

    def num(s):
        return int(s.split('.')[0].split('_')[-1])

    global data
    files = os.listdir(anim_folder)
    files.sort(key=num)
    print(files)
    data = []
    for f in files: 
        data.append(np.transpose(pd.read_csv(anim_folder + f, header=None).values))
    

load_files()

fig, axes = plt.subplots(nrows=1, figsize=(10, 4))

# image = axes.imshow(all_data[0], vmin=0, vmax=1)
vmin = 1e100
vmax = 0

for d in data:
    if vmax < d.max().max(): vmax = d.max()
    if vmin > d.min().min(): vmin = d.min() 

image = axes.imshow(data[0], cmap='plasma', origin='lower', vmin=vmin, vmax=vmax)

divider0 = make_axes_locatable(axes)
cax0 = divider0.append_axes("right", size="1%", pad=0.1)
fig.colorbar(image, ax=axes, cax=cax0)


animation = FuncAnimation(
    # Your Matplotlib Figure object
    fig,
    # The function that does the updating of the Figure
    animate,
    # Frame information (here just frame number)
    np.arange(25,50),
    # Extra arguments to the animate function
    fargs=[],
    # Frame-time in ms; i.e. for a given frame-rate x, 1000/x
    interval=1000 / 10
)

animation.save("anim.mp4", dpi=300)

# if __name__=='__main__':
    


