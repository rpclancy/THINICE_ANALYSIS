import numpy as np
import sys

from IPython.display import clear_output

def regrid_spherical(info, all_systems, var_anom, var_dates,
                     latmesh, lonmesh, ACorTPV, var_source):
    """
    For each AC/TPV
        Get distance and bearing to every ERA5 grid point and the variable
        Sometimes only those where sea ice is present and in hemis)
        Composit onto grid, keep composit and note which grid cells have data.
    """
    
    # Unpack input class
    (nc_file, var_name, var_date_name, t_ind, scale, months, amp_perc, max_lat, max_d, step_size_a,
     step_size_d, SIC_only, nonanom, regions, change, dstart, doff, ystart, yend) = unpack_spatial_info(info)
    
    # Get AC or TPV data
    if ACorTPV == 'TPV':
        sys_year = all_systems[:,0]
        sys_month = all_systems[:,1]
        sys_day = all_systems[:,2]
        sys_hour = all_systems[:,3]
        sys_latitude = all_systems[:,4]
        sys_longitude = all_systems[:,5]
    elif ACorTPV == 'AC':
        sys_year = all_systems[:,2]
        sys_month = all_systems[:,3]
        sys_day = all_systems[:,4]
        sys_hour = all_systems[:,5]
        sys_latitude = all_systems[:,13]
        sys_longitude = all_systems[:,14]
    if ACorTPV == 'SCAC':
        sys_year = all_systems[:,0]
        sys_month = all_systems[:,1]
        sys_day = all_systems[:,2]
        sys_hour = all_systems[:,3]
        sys_latitude = all_systems[:,4]
        sys_longitude = all_systems[:,5]
    
    # Get dates/time
    var_year = var_dates[0,:]
    var_month = var_dates[1,:]
    var_day = var_dates[2,:]
    if var_source == 'ERA5':
        var_hour = var_dates[3,:]
    
    # Create interp grid
    a_il = np.arange(0, 361 - step_size_a, step_size_a)
    a_iu = a_il + step_size_a
    d_il = np.arange(0, max_d+1 - step_size_d, step_size_d)
    d_iu = d_il + step_size_d 
    
    # Set up variables to collect spherical regriddings
    var_grid_all = []
    grid_count_all = []
    
    # For each AC or TPV
    max_count = len(sys_year)
    for i in range(0, len(sys_year)):
        print(i, 'out of', max_count)
        
        # Skip years out of range of analysis
        if sys_year[i] < np.min(var_year):
            continue
        
        # Get lat/lon of individual system
        sys_latitude_sub = sys_latitude[i]
        sys_longitude_sub = sys_longitude[i]
        
        # Find index of matching variable dates
        if var_source == 'ERA5':
            var_ind = (var_year==sys_year[i]) & (var_month==sys_month[i]) & (var_day==sys_day[i]) & (var_hour==sys_hour[i])
        else:
            var_ind = (var_year==sys_year[i]) & (var_month==sys_month[i]) & (var_day==sys_day[i])
        
        if sum(var_ind) < 1:
            if var_source == 'ERA5':
                print('No data:', sys_year[i], sys_month[i], sys_day[i], sys_hour[i])
            else:
                print('No data:', sys_year[i], sys_month[i], sys_day[i])
            continue
        
        # Getting var_anom for that date
        var_anom_sub = np.squeeze(var_anom[:,:,var_ind])
        
        # Getting distances in km and bearing to every grid cell
        distances, bearing = rhumb_np(sys_longitude_sub, sys_latitude_sub, lonmesh, latmesh)
        
        # Setting to nan points in different hemisphere (meridionally)
        in_hemis = find_in_hemisphere(sys_longitude_sub, lonmesh)
        var_anom_sub[in_hemis!=1] = np.nan
        
        # Setting to nan points outside of max_d of system
        var_anom_sub[distances>=max_d] = np.nan
        
        # Flatten
        distances = np.ndarray.flatten(distances)
        bearing = np.ndarray.flatten(bearing)
        var_anom_sub = np.ndarray.flatten(var_anom_sub)
        
        # Select only non-nan values
        distances = distances[np.isfinite(var_anom_sub)]
        bearing = bearing[np.isfinite(var_anom_sub)]
        var_anom_sub = var_anom_sub[np.isfinite(var_anom_sub)]
        
        var_grid = np.zeros([a_il.size, d_il.size])
        grid_count = np.zeros([a_il.size, d_il.size])
        
        # For each distance/bearing combo average values together
        for i in range(0, a_il.size):
            for j in range(0, d_il.size):
                cond_1 = bearing>a_il[i]
                cond_2 = bearing<a_iu[i]
                cond_3 = distances>d_il[j]
                cond_4 = distances<d_iu[j]
        
                grid_ind = (np.column_stack((cond_1, cond_2, cond_3, cond_4)).all(axis=1))
                grid_count[i,j] = np.count_nonzero(grid_ind)
                var_grid_sub = var_anom_sub[grid_ind]
                var_grid[i,j] = np.nanmean(var_grid_sub, dtype=np.float64) 
                #dtype=np.float64 is important for accuracy if using float32 or float16 as inputs
        
        # Add to list of regridded variables.
        var_grid_all.append(var_grid)
        grid_count_all.append(grid_count)
        
    return var_grid_all, grid_count_all, a_il, a_iu, d_il, d_iu
    

def unpack_spatial_info(info):
    """ Unpack input class """
    nc_file = info.nc_file
    var_name = info.var_name
    var_date_name = info.var_date_name
    t_ind = info.t_ind
    scale = info.scale
    months = info.months
    amp_perc = info.amp_perc
    max_lat = info.max_lat
    max_d = info.max_d
    step_size_a = info.step_size_a
    step_size_d = info.step_size_d
    SIC_only = info.SIC_only
    nonanom = info.nonanom
    regions = info.regions
    change = info.change
    dstart = info.dstart
    doff = info.doff
    ystart = info.ystart
    yend = info.yend
    return (nc_file, var_name, var_date_name, t_ind, scale, months,
            amp_perc, max_lat, max_d, step_size_a, step_size_d, SIC_only,
            nonanom, regions, change, dstart, doff, ystart, yend)

def rhumb_np(lon1, lat1, lons, lats):
    """
    Defining distance and bearing formula for rhumb line
    (based on https://www.movable-type.co.uk/scripts/latlong.html)
    """
    lon1, lat1, lons, lats = map(np.radians, [lon1, lat1, lons, lats])
    
    dlats = lats - lat1
    dlons = lons - lon1
    dlons[dlons < -np.pi] = dlons[dlons < -np.pi] + 2*np.pi
    dlons[dlons > np.pi] = dlons[dlons > np.pi] - 2*np.pi
    
    #Δψ = ln( tan(π/4 + φ2/2) / tan(π/4 + φ1/2) )
    dpsi = np.log( np.tan(np.pi/4 + lats/2) / np.tan(np.pi/4 + lat1/2) )
    
    #q = Δφ/Δψ
    q=np.zeros_like(dlats)
    q[dlats!=0] = dlats[dlats!=0]/dpsi[dlats!=0]
    #or q=cosφ for E-W line
    q[dlats==0] = np.cos(lat1)
    
    #Calculate distance
    #d = √(Δφ² + q²⋅Δλ²) ⋅ R
    #Δλ is taking shortest route (<180°)
    R = 6367 #roughly radius of earth in km
    d = np.sqrt(dlats*dlats + q*q*dlons*dlons) * R 
    
    #Calculate bearing
    #θ = atan2(Δλ, Δψ)
    bearing = np.arctan2(dlons, dpsi)
    bearing = np.degrees(bearing)
    bearing = (bearing + 360) %360
    
    return d, bearing

def find_in_hemisphere(system_lon, lonmesh):
    """
    Finds which ecmwf grid cells are the wrong side of the pole
    Bearing gets kinda meaningless in that case (e.g. can be south but point N)
    Gives 1 for cells right side of pole, 0 for 
    """
    A = system_lon
    B = lonmesh
    in_hemis = (abs(B-A)<90)*1 + (abs(B-A+360)<90)*1 + (abs(B-A-360)<90)*1
    in_hemis = np.sign(in_hemis) #might be redundant
    
    return(in_hemis)
