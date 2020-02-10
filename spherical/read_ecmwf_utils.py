import numpy as np

from netCDF4 import Dataset

def read_nc_var(nc_file, var_name, var_date_name):
    """ 
    Reads in data and associated dates from netcdf file.
    N.B. These netcdf files were created through matlab workflow.
    """
    
    fh = Dataset(nc_file, mode='r')
    var = fh.variables[var_name][:]
    var_dates = fh.variables[var_date_name][:]
    
    print(var_name,': ', var.shape)
    print(var_date_name,': ', var_dates.shape)
    
    return var, var_dates;

def get_y_mo_d(var_dates):
    """ Gets year, month, day vectors from dates extracted using read_nc_var"""
    v_year = var_dates[0,:]
    v_month = var_dates[1,:]
    v_day = var_dates[2,:]
    
    return v_year, v_month, v_day;

def get_ecmwf_lat_lon(nc_file):
    """ Gets ecmwf lats and lons from netcdf file and creates meshgrids"""
    from netCDF4 import Dataset
    
    fh = Dataset(nc_file, mode='r')

    latitude_ecmwf = fh.variables['latitude_ecmwf'][:]
    longitude_ecmwf = fh.variables['longitude_ecmwf'][:]

    lonmesh_ecmwf,latmesh_ecmwf = np.meshgrid(longitude_ecmwf,latitude_ecmwf)

    print('latitude_ecmwf: ', latitude_ecmwf.shape)
    print('longitude_ecmwf: ', longitude_ecmwf.shape)
    
    return latitude_ecmwf, longitude_ecmwf, latmesh_ecmwf, lonmesh_ecmwf;

def calculate_clim_anoms(var, var_dates):
    """
    Calculates climatology of a (lat x lon x time) variable and anomalies from this.
    Also calculates anomaly scaled by climatology standard deviation at that grid point:
        i.e. how unusual an anomaly is.
    """
    d_counts=[]
    var_clim = np.zeros_like(var)
    var_climstd = np.zeros_like(var)
    for m in range(1,13): #for each month
        mo_ind = (var_dates[1,:]==m)
        day_options = np.unique(var_dates[2,mo_ind])
    
        #print(day_options) #for diagnostics 
        for d in range(0,np.size(day_options)): #for each possible day
            d_ind = (mo_ind) & (var_dates[2,:]==day_options[d])

            var_days = var[:,:,d_ind]
            var_daysav = np.nanmean(var_days,2)
            var_daysstd = np.nanstd(var_days,2)
        
            var_clim[:,:,d_ind] = np.transpose(np.tile(var_daysav,(np.sum(d_ind),1,1)),(1,2,0))
            var_climstd[:,:,d_ind] = np.transpose(np.tile(var_daysstd,(np.sum(d_ind),1,1)),(1,2,0))
        
            d_counts.append(np.sum(d_ind)) #this is just for diagnostics
        
    var_anom = var - var_clim
    var_anom_scaled = var_anom/var_climstd
    
    return var_anom, var_anom_scaled;

def calculate_climatology(var, var_dates):
    """
    Calculates climatology of a (lat x lon x time) variable.
    n.b. untested as of this moment. not sure if I'll use it in the end even.
    """
    d_counts=[]
    var_clim = np.zeros_like(var)
    for m in range(1,13): #for each month
        mo_ind = (var_dates[1,:]==m)
        day_options = np.unique(var_dates[2,mo_ind])
    
        #print(day_options) #for diagnostics 
        for d in range(0,np.size(day_options)): #for each possible day
            d_ind = (mo_ind) & (var_dates[2,:]==day_options[d])

            var_days = var[:,:,d_ind]
            var_daysav = np.nanmean(var_days,2)
        
            var_clim[:,:,d_ind] = np.transpose(np.tile(var_daysav,(np.sum(d_ind),1,1)),(1,2,0))
        
            d_counts.append(np.sum(d_ind)) #this is just for diagnostics
    
    return var_clim;

def generate_sic_monthly():
    """
    Finds grid cells where sea ice is present at any point in each month
    """
    # Set variable
    nc_file = '/home/disk/sipn/rclancy/ecmwf/pf/predictability/SIC/SIC.nc'
    var_name = 'SIC'
    var_date_name = 'SIC_dates'

    # Get variables
    var, var_dates =read_nc_var(nc_file, var_name, var_date_name)
    # Get associated dates
    var_year, var_month, var_day = get_y_mo_d(var_dates)

    SIC_present=np.zeros([12, 20, 240])
    for m in range(1,13):
        mo_ind = np.array(var_month==m)
        SIC_present[m-1,:,:] = np.nanmean(var[:,:,mo_ind],2)

    SIC_present[SIC_present>0]=1
    return SIC_present;