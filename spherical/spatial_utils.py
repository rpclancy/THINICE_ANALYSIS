import numpy as np

from IPython.display import clear_output

import read_ecmwf_utils as reu
import read_TPV_utils as rtu
import plot_utils as pu


def haversine_np(lon1, lat1, lons, lats):
    """ Defining distance formula (assumes spherical earth, probably fine?) """
    lon1, lat1, lons, lats = map(np.radians, [lon1, lat1, lons, lats])

    dlons = lons - lon1
    dlats = lats - lat1

    a = np.sin(dlats/2.0)**2 + np.cos(lat1) * np.cos(lats) * np.sin(dlons/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km

def azimuth_np(lon1, lat1, lons, lats):
    """ Defining azimuth formula (based on https://www.movable-type.co.uk/scripts/latlong.html) """
    lon1, lat1, lons, lats = map(np.radians, [lon1, lat1, lons, lats])
    
    dlons = lons - lon1
    dlats = lats - lat1
    
    bvar_1 = np.sin(dlons)*np.cos(lats)
    bvar_2 = np.cos(lat1)*np.sin(lats) - np.sin(lat1)*np.cos(lats)*np.cos(dlons)
    bearing = np.arctan2(bvar_1,bvar_2)
    
    bearing = np.degrees(bearing)
    bearing = (bearing + 360) %360
    return bearing

def rhumb_np(lon1, lat1, lons, lats):
    """
    Defining distance and bearing formula for rhumb line
    (based on https://www.movable-type.co.uk/scripts/latlong.html)
    
    n.b. haven't considered here the implications of going across north pole, but could maybe
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

def find_in_hemisphere(TPV_lon, lonmesh_ecmwf):
    """
    Finds which ecmwf grid cells are the wrong side of the pole
    Bearing gets kinda meaningless in that case (e.g. can be south but point N)
    Gives 1 for cells right side of pole, 0 for 
    """
    A = TPV_lon
    B = lonmesh_ecmwf 
    in_hemis = (abs(B-A)<90)*1 + (abs(B-A+360)<90)*1 + (abs(B-A-360)<90)*1
    in_hemis = np.sign(in_hemis) #might be redundant
    
    return(in_hemis)

def var_TPV_spatial(all_TPVs_A, var_dates, var_anom, latmesh_ecmwf, lonmesh_ecmwf):
    """
    Get variable anomaly for each grid point and E/W and N/S distance to TPV for that grid point
    """
    # set up variables
    # 20 by 240 is to match sizes of latmesh_ecmwf and lonmesh_ecmwf
    c = 0 #tpv count variable
    d = np.zeros([20,240,15000])
    b = np.zeros([20,240,15000])
    var_spatial = np.zeros([20,240,15000])
    lats_spatial = np.zeros([20,240,15000])
    amp_spatial = np.zeros([20,240,15000])
    in_hemis = np.zeros([20,240,15000])
    SIC_present_spatial = np.zeros([20,240,15000])
    
    var_year = var_dates[0,:]
    var_month = var_dates[1,:]
    var_day = var_dates[2,:]
    
    TPV_year = all_TPVs_A[:,0]
    TPV_month = all_TPVs_A[:,1]
    TPV_day = all_TPVs_A[:,2]
    TPV_latitude = all_TPVs_A[:,4]
    TPV_longitude = all_TPVs_A[:,5]
    TPV_amplitude = all_TPVs_A[:,7]
    
    # Get where sea ice was present in each month
    SIC_present = reu.generate_sic_monthly()
    
    # For each variable date
    for i in range(0,np.size(var_dates,1)):
        print(i, 'out of', np.size(var_dates,1))
        # See if there are any TPVs
        TPV_inds = (TPV_year==var_year[i]) & (TPV_month==var_month[i]) & (TPV_day==var_day[i])
        TPV_count = np.sum(TPV_inds)
        if TPV_count>0:
            #getting lats and lons of these TPVs
            TPV_latitude_sub = TPV_latitude[TPV_inds]
            TPV_longitude_sub = TPV_longitude[TPV_inds]
            TPV_amplitude_sub = TPV_amplitude[TPV_inds]
            #getting var_anom for that date
            var_anom_sub = var_anom[:,:,i]
            
            #getting where ice was present that month
            SIC_present_sub = SIC_present[int(var_month[i]-1),:,:]
            
            #getting distances in km to every grid cell from each TPV
            for j in range(0,np.size(TPV_longitude_sub,0)):
                distances, bearing = rhumb_np(TPV_longitude_sub[j], TPV_latitude_sub[j],
                                         lonmesh_ecmwf, latmesh_ecmwf)
                
                d[:,:,c] = distances
                b[:,:,c] = bearing
                var_spatial[:,:,c] = var_anom_sub
                lats_spatial[:,:,c] = latmesh_ecmwf
                amp_spatial[:,:,c] = TPV_amplitude_sub[j]
                SIC_present_spatial[:,:,c] = SIC_present_sub
                
                in_hemis[:,:,c] = find_in_hemisphere(TPV_longitude_sub[j], lonmesh_ecmwf)
                
                c+=1
    
    #trimming trailing zeros
    d = d[:,:,0:c]
    b = b[:,:,0:c]
    var_spatial = var_spatial[:,:,0:c]
    lats_spatial = lats_spatial[:,:,0:c]
    amp_spatial = amp_spatial[:,:,0:c]
    in_hemis = in_hemis[:,:,0:c]
    SIC_present_spatial = SIC_present_spatial[:,:,0:c]

    return d, b, var_spatial, lats_spatial, amp_spatial, in_hemis, SIC_present_spatial;

def select_lats(d, b, var_spatial, lats_spatial, amp_spatial, in_hemis, max_lat = 90):
    """ Selects values south of a given latitude as E/W is poorly defined near 90N """
    d = d[lats_spatial<=max_lat]
    b = b[lats_spatial<=max_lat]
    var_spatial = var_spatial[lats_spatial<=max_lat]
    amp_spatial = amp_spatial[lats_spatial<=max_lat]
    in_hemis = in_hemis[lats_spatial<=max_lat]
    lats_spatial = lats_spatial[lats_spatial<=max_lat] #make sure is last
    
    return d, b, var_spatial, lats_spatial, amp_spatial, in_hemis;

def select_amps(d, b, var_spatial, lats_spatial, amp_spatial, in_hemis, amp_perc = 100):
    """ Selects values related to the highest percentile (amp_perc) of TPV amplitudes """
    #Get threshold value
    amp_co = np.percentile(amp_spatial, 100-amp_perc)
    d = d[amp_spatial>=amp_co]
    b = b[amp_spatial>=amp_co]
    var_spatial = var_spatial[amp_spatial>=amp_co]
    lats_spatial = lats_spatial[amp_spatial>=amp_co]
    in_hemis = in_hemis[amp_spatial>=amp_co]
    amp_spatial = amp_spatial[amp_spatial>=amp_co] #make sure is last
    
    return d, b, var_spatial, lats_spatial, amp_spatial, in_hemis;

def select_SIC_only(d, b, var_spatial, lats_spatial, amp_spatial, in_hemis, SIC_present_spatial):
    """Selects only values over where there was sea ice at any time in that month"""
    d = d[SIC_present_spatial==1]
    b = b[SIC_present_spatial==1]
    var_spatial = var_spatial[SIC_present_spatial==1]
    lats_spatial = lats_spatial[SIC_present_spatial==1]
    in_hemis = in_hemis[SIC_present_spatial==1]
    amp_spatial = amp_spatial[SIC_present_spatial==1]
    
    return d, b, var_spatial, lats_spatial, amp_spatial, in_hemis;
                                            

def grid_TPV_spatial(d, b, var_spatial, lats_spatial, amp_spatial, in_hemis, SIC_present_spatial,
                     amp_perc=100, max_lat=90, max_d=1000, step_size_a=30, step_size_d=100, SIC_only=False):
    """
    Copositing variable anomalies on grid of distances from TPV
    #max_d = maximum distance from tpv in km
    #SIC_only = True/False for if only do anlaysis over sea ice
    """
    
    # select only points where there was sea ice that month if necessary
    if SIC_only:
        d, b, var_spatial, lats_spatial, amp_spatial, in_hemis = select_SIC_only(d, b, var_spatial,
                                            lats_spatial, amp_spatial, in_hemis, SIC_present_spatial)
        
    # select only points below prescribed max latitude if necessary
    if max_lat != 90:
        d, b, var_spatial, lats_spatial, amp_spatial, in_hemis = select_lats(d, b, var_spatial,
                                                    lats_spatial, amp_spatial, in_hemis, max_lat)
        
    # select upper fraction of TPV amplitudes if necessary
    if amp_perc !=100:
        d, b, var_spatial, lats_spatial, amp_spatial, in_hemis = select_amps(d, b, var_spatial,
                                                    lats_spatial, amp_spatial, in_hemis, amp_perc)
    
    d_sub = np.ndarray.flatten(d)
    b_sub = np.ndarray.flatten(b)
    var_sub = np.ndarray.flatten(var_spatial)
    in_hemis = np.ndarray.flatten(in_hemis)
    
    # select only points on right side of pole
    b_sub = b_sub[in_hemis==1]
    d_sub = d_sub[in_hemis==1]
    var_sub = var_sub[in_hemis==1]

    # select only points within max_d of tpv
    b_sub = b_sub[d_sub<=max_d]
    var_sub = var_sub[d_sub<=max_d]
    d_sub = d_sub[d_sub<=max_d]
    
    # create interp grid
    a_il = np.arange(0, 361 - step_size_a, step_size_a)
    a_iu = a_il + step_size_a
    d_il = np.arange(0, max_d+1 - step_size_d, step_size_d)
    d_iu = d_il + step_size_d

    var_grid = np.zeros([a_il.size, d_il.size])
    grid_count = np.zeros([a_il.size, d_il.size])

    for i in range(0, a_il.size):
        clear_output(wait=True)
        print(i, ' out of', a_il.size)
        for j in range(0, d_il.size):
            cond_1 = b_sub>a_il[i]
            cond_2 = b_sub<a_iu[i]
            cond_3 = d_sub>d_il[j]
            cond_4 = d_sub<d_iu[j]
            #cond_5 = ~np.isnan(var_sub) #don't think is needed.
        
            grid_ind = (np.column_stack((cond_1, cond_2, cond_3, cond_4)).all(axis=1))
            grid_count[i,j] = np.count_nonzero(grid_ind)
            var_grid_sub = var_sub[grid_ind]
            var_grid[i,j] = np.nanmean(var_grid_sub)
            
    return var_grid, grid_count, a_il, a_iu, d_il, d_iu;

def var_TPV_full_spatial(info):
    """
    Same as var_TPV_full_spatial, but trying to rework so that it takes a class instance as input
    """
    # Unpack input class
    nc_file = info.nc_file
    var_name = info.var_name
    var_date_name = info.var_date_name
    t_ind = info.t_ind
    scale = info.scale
    amp_perc = info.amp_perc
    max_lat = info.max_lat
    max_d = info.max_d
    step_size_a = info.step_size_a
    step_size_d = info.step_size_d
    SIC_only = info.SIC_only
    nonanom = info.nonanom
    
    # Read data
    var, var_dates = reu.read_nc_var(nc_file, var_name, var_date_name)

    # Get associated dates
    var_year, var_month, var_day = reu.get_y_mo_d(var_dates)
    
    # Get lat/lon coords
    latitude_ecmwf, longitude_ecmwf, latmesh_ecmwf, lonmesh_ecmwf = reu.get_ecmwf_lat_lon(nc_file)
    
    # For ensemble spread variables: get specified timestep
    if t_ind is not False:
        var = var[t_ind,:,:,:]
    
    # Calculate variable anomalies
    var_anom, var_anom_scaled = reu.calculate_clim_anoms(var, var_dates)
    
    # Load TPV data
    all_TPVs_A = np.loadtxt('../all_TPVs_A_matched')
    colnames_TPVs = open('../colnames_TPVs.ascii').read().split()
    
    # Select just JJA TPVs north of 65N from 1979 to 2016
    all_TPVs_A = rtu.select_JJA_65N_1979to2016(all_TPVs_A)
    
    # Get variable anomaly for each grid point and
    # azimuth/distance to TPV for that grid point (using rhumb lines)
    if scale:
        var_anom = var_anom_scaled
    elif nonanom:
        var_anom = var
        
    d, b, var_anom_spatial, lats_spatial, amp_spatial, in_hemis, SIC_present_spatial= var_TPV_spatial(
                                        all_TPVs_A, var_dates, var_anom, latmesh_ecmwf, lonmesh_ecmwf)
    
    # Copositing variable anomalies on grid of distances from TPV
    var_grid, grid_count, a_il, a_iu, d_il, d_iu = grid_TPV_spatial(d, b, var_anom_spatial, lats_spatial,
            amp_spatial, in_hemis, SIC_present_spatial, amp_perc, max_lat, max_d, step_size_a, step_size_d, SIC_only);

    # Getting grid midpoints
    a_im = (a_il + a_iu)/2
    d_im = (d_il + d_iu)/2
    
    return var_grid, grid_count, a_im, d_im;