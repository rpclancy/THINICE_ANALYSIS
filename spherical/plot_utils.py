import matplotlib.pyplot as plt
import numpy as np

def spatial_plot_format(var_name, scale):
    """ Formats TPV spatial plots consistently """
    title_str = "TPV effect on " + var_name
    plt.title(title_str)
    
    cbar_str = var_name
    if scale:
        cbar_str = cbar_str + " (stdv)"
    plt.colorbar(fraction=0.023, pad=0.04, label = cbar_str)
    
    if var_name != 'grid_count':
        A, B = plt.gci().get_clim()
        cmax=np.abs(np.array([A, B])).max()
        plt.clim(-cmax, cmax)
    else:
        A, B = plt.gci().get_clim()
        cmax=np.abs(np.array([A, B])).max()
        plt.clim(0, cmax/3)
        plt.title('Grid box count')
    
    plt.axhline(y=0, color='k', linestyle='-')
    plt.axvline(x=0, color='k', linestyle='-')
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    
def contourf_tpv_spatial(var_grid, x_im, y_im, var_name, scale):
    """ Plotting variable relative to TPV position"""
    fig=plt.figure(figsize=(6, 5), dpi= 80, facecolor='w', edgecolor='k')
    plt.contourf(x_im, y_im, np.transpose(var_grid), cmap='seismic')
    spatial_plot_format(var_name, scale)
    
def pcolor_tpv_spatial(var_grid, x_im, y_im, var_name, scale):
    """ Plotting variable relative to TPV position"""
    fig=plt.figure(figsize=(6, 5), dpi= 80, facecolor='w', edgecolor='k')
    plt.pcolor(x_im, y_im, np.transpose(var_grid), cmap='seismic')
    spatial_plot_format(var_name, scale)  

def pcolor_gridcount_spatial(grid_count, x_im, y_im):    
    """ Plotting number of data points in each grid square"""
    fig=plt.figure(figsize=(6, 5), dpi= 80, facecolor='w', edgecolor='k')
    plt.pcolor(x_im, y_im, np.transpose(grid_count), cmap='jet')
    spatial_plot_format('grid_count', scale=False)
