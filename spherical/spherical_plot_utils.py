import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def spherical_plot_format(var_name, scale):
    """ Formats TPV spherical plots consistently """
    title_str = "TPV effect on " + var_name
    plt.title(title_str)

    cbar_str = var_name
    if scale:
        cbar_str = cbar_str + " (stdv)"
    plt.colorbar(fraction=0.023, pad=0.04, label = cbar_str)
    
    A, B = plt.gci().get_clim()
    cmax=np.abs(np.array([A, B])).max()
    
    if var_name != 'grid_count':
        plt.clim(-cmax, cmax)
    else:
        plt.clim(0, cmax)
        plt.title('Grid box count')

def pcolor_tpv_spherical(var_grid, step_size_d, max_d, step_size_a, var_name, scale):
    """ Plotting variable relative to TPV position"""
    r, th = np.meshgrid(np.arange(0,max_d+1,step_size_d), np.radians(np.arange(0,361,step_size_a)))    
    fig=plt.figure(figsize=(6, 5), dpi= 80, facecolor='w', edgecolor='k')
    ax = plt.subplot(projection="polar")
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    plt.pcolormesh(th, r, var_grid, cmap='seismic')
    spherical_plot_format(var_name, scale)  

def pcolor_gridcount_spherical(grid_count, step_size_d, max_d, step_size_a): 
    """ Plotting number of data points in each grid square"""
    r, th = np.meshgrid(np.arange(0,max_d+1,step_size_d), np.radians(np.arange(0,361,step_size_a)))    
    fig=plt.figure(figsize=(6, 5), dpi= 80, facecolor='w', edgecolor='k')
    ax = plt.subplot(projection="polar")
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    plt.pcolormesh(th, r, grid_count, cmap='jet')
    spherical_plot_format('grid_count', scale=False)