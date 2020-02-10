import numpy as np

def save_spatial_np(info, var_grid, grid_count, x_im, y_im):
    """
    Saves results of spatial analysis in numpy arrays for later use 
    Probably superceded by save_spherical_np
    """    
    #Unpack inputs
    var_name = info.var_name
    t_ind = info.t_ind
    scale = info.scale
    max_lat = info.max_lat
    max_d = info.max_d
    step_size = info.step_size
    day_ind = info.day_ind
    
    d = False
    if t_ind is not False:
        d = day_ind[t_ind]
    
    file_name = f"{var_name}_scaled={scale}_forecastd={d}_maxlat={max_lat}_maxd={max_d}_step={step_size}"
    file_name = 'result_data/' + file_name + '.npz'
    
    if var_name == 'PSL':
        np.savez(file_name, PSL_grid=var_grid, grid_count=grid_count, x_im=x_im, y_im=y_im)
    elif var_name == 'Z500':
        #We'll fill this in later
        pass

def save_spherical_np(info, var_grid, grid_count):
    """ Saves results of spherical analysis in numpy arrays for later use """ 
    #Unpack inputs
    var_name = info.var_name
    t_ind = info.t_ind
    scale = info.scale
    amp_perc = info.amp_perc
    max_lat = info.max_lat
    max_d = info.max_d
    step_size_d = info.step_size_d
    step_size_a = info.step_size_a
    day_ind = info.day_ind

    d = False
    if t_ind is not False:
        d = day_ind[t_ind]
        
    file_name = f"{var_name}_scaled={scale}_forecastd={d}_ampperc={amp_perc}_maxlat={max_lat}_maxd={max_d}_stepd={step_size_d}_stepa={step_size_a}"
    file_name = 'result_data/' + file_name + '.npz'
    
    if var_name == 'PSL' or var_name == 'PSL_RMSE' or var_name == 'PSL_stdv':
        np.savez(file_name, PSL_grid=var_grid, grid_count=grid_count,
                 step_size_a=step_size_a, step_size_d=step_size_d, max_d=max_d)
    elif var_name == 'Z500' or var_name == 'Z500_RMSE' or var_name == 'Z500_stdv':
        np.savez(file_name, Z500_grid=var_grid, grid_count=grid_count,
                 step_size_a=step_size_a, step_size_d=step_size_d, max_d=max_d)
    elif var_name == 'SIC' or var_name == 'SIC_RMSE' or var_name == 'SIC_stdv':
        np.savez(file_name, SIC_grid=var_grid, grid_count=grid_count,
                 step_size_a=step_size_a, step_size_d=step_size_d, max_d=max_d)
    elif var_name == 'U' or var_name == 'U_RMSE' or var_name == 'U_stdv':
        np.savez(file_name, U_grid=var_grid, grid_count=grid_count,
                 step_size_a=step_size_a, step_size_d=step_size_d, max_d=max_d)
    elif var_name == 'V' or var_name == 'V_RMSE' or var_name == 'V_stdv':
        np.savez(file_name, V_grid=var_grid, grid_count=grid_count,
                 step_size_a=step_size_a, step_size_d=step_size_d, max_d=max_d)
    elif var_name == 'something_else':
        #We'll fill this in later
        pass