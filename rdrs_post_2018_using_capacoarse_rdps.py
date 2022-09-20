#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Description

"""

__author__      = 'Daren Tsang/Dikra Khedhaouiria/Milena Dimitrijevic'
__copyright__   = 'Copyright 2022, RDPA format as RDRS'
__license__     = '{license}'
__version__     = '0.1.0'
__maintainer__  = '{Dikra Khedhaouiria}'
__email__       = '{dikraa.khedhaouiria@ec.gc.ca}'
__status__      = 'Dev'


import netCDF4 as nc4
import numpy as np    
import time        
from configparser import ConfigParser
import math
from os.path import join, exists
import os
from datetime import date, timedelta
import sys
import logging
from scipy.spatial import KDTree
import warnings; warnings.filterwarnings(action='ignore')
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

def read_data(RDPS_directory, CaPA_directory, RDRS_template, start_date):
    """
    Reads in the input netCDF files
    Parameters
    ----------
    RDPS_directory : str
        path of the rdps files.
    CaPA_directory : str
        path of the capa files.
    RDRS_template : str
        name of the rdrs file used as template.
    start_date : str
        start date in the format of YYYYMMDD

    Returns
    -------
    RDPS_files : list
        list of the rdps files read as netcdf.
    CaPA_files : list
        list of the capa files read as netcdf..
    RDRS_file : netcdf
        netcdf object with the information of the rdrs file used as template.
    namefile: dict
        
    """
 
    # turn the start_date datetime object into a string
    start_date_string = start_date.strftime("%Y%m%d")
    
    # also get the string for the date following the start date
    next_date         = start_date + timedelta(days = 1)
    next_date_string  = next_date.strftime("%Y%m%d")
    
    # create empty list to store the RDPS, CaPA and the name of files
    RDPS_files  = []
    namefile    = {}
    CaPA_files  = []
    
    # list of needed files from RDPS and CaPA
    namefile['RDPS'] = [join(RDPS_directory, start_date_string + "06.nc"),
                        join(RDPS_directory, start_date_string + "12.nc"),
                        join(RDPS_directory, start_date_string + "18.nc"),
                        join(RDPS_directory, next_date_string  + "00.nc")]
    
    namefile['CaPA'] = [join(CaPA_directory, start_date_string + "18.nc"),
                        join(CaPA_directory, next_date_string + "00.nc"),
                        join(CaPA_directory, next_date_string + "06.nc"),
                        join(CaPA_directory, next_date_string + "12.nc")]
            
    # Check if all RDPS files are present
    for file in namefile['RDPS']: 
        try:
            RDPS_files.append(nc4.Dataset(file, 'r'))
            logging.info(f"RDPS file {file} is present --> OK")
        except:
            sys.stderr.write(f"PROBLEM opening (or non existing) RDPS file: {file}")
        
    
    # Check if all CaPA files are present
    for file in namefile['CaPA']: 
        try:
            CaPA_files.append(nc4.Dataset(file, 'r'))    
            logging.info(f"CaPA file {file} is present --> OK")
        except:
            sys.stderr.write(f"PROBLEM opening (or non existing) CaPA file: {file}")
        
    # read the RDRS file to use as a template
    RDRS_file = nc4.Dataset(RDRS_template, 'r')
    
    #print (print(RDRS_file.__dict__))
    #name_product_RDRS = str(RDRS_file.__dict__['product'])

    return RDPS_files, CaPA_files, RDRS_file, namefile

    
def check_vars(RDPS_files, var_list, nameprodfile):
    """
    Checks to make sure all variables in the list are present in given RDPS files

    Parameters
    ----------
    RDPS_files : list
        DESCRIPTION.
    var_list : list
        list of the variable name.
    nameprodfile : dict
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    compt = 0
    for file in RDPS_files: 
        name_file     = nameprodfile['RDPS'][compt]
        RDPS_var_list = list(file.variables.keys())
        
        lst_var_bl    = [x in RDPS_var_list for x in var_list]
        
        # checks if var_list is a subset of RDPS_variables list
        if(all(lst_var_bl)):
            pass
        else:
            
            missingvars = list(np.where(~np.array(lst_var_bl))[0])
            var_list_np = np.array(var_list)
            
            sys.exit(f"Variables {var_list_np[missingvars]} are missing in {name_file}")
            
        compt = compt + 1
    logging.info(f"All variables {var_list} are present in each RDPS files --> OK")
    
 
def check_grid(prod_files, nameproddict, nameprod):
    """
    Checks to make sure all RDPS_files have the same grid definition, and all CaPA files have the same grid definition

    Parameters
    ----------
    prod_files : list 
        list of netcdf4 datasets (RDPS or CaPA).
    nameproddict : dict
        dictionnary with the corresponding file name.
    nameprod : str
        name of the product (RDPS or CaPA).

    Returns
    -------
    None.

    """
    
    base_PROD_grid, base_PROD_lat, base_PROD_lon = [], [], []
    for file in prod_files:
        base_PROD_grid.append(file.variables["lat"].eccc_grid_definition)
        base_PROD_lat.append(file.variables["lat"][:])
        base_PROD_lon.append(file.variables["lon"][:])    
    
    # check the grid definition   
    if len(np.unique(base_PROD_grid)) == 1:
        pass
    else:
        sys.exit(f"{nameprod} files must be on the same grid: {nameproddict[nameprod]}")
    
    logging.info(f"All {nameprod} files have the same grid")   
    
    # check the domain definition
    for i in range(len(base_PROD_grid)):
        for j in range(i,len(base_PROD_grid)):
            
            equal_or_not_PROD_lat = np.ma.allequal(base_PROD_lat[i], base_PROD_lat[j])
            equal_or_not_PROD_lon = np.ma.allequal(base_PROD_lon[i], base_PROD_lon[j])
            if equal_or_not_PROD_lat and equal_or_not_PROD_lon:
                pass
            else:
                sys.exit(f"{nameprod} file: {nameproddict[nameprod]} has a different domain")
                
    logging.info(f"All {nameprod} files have the same domain")   
            

def check_time(RDPS_files, required_time, nameprodfile):
    """
    Checks to make sure the RDPS files have the required 7 hours of data, 
    and that each file contains the same forecast horizons
    
    Parameters
    ----------
    RDPS_files : list
         list of netcdf4 datasets (RDPS).
    required_time : list
        list of the required time steps (lead time 6-12)
    nameprodfile : dict
        information on the name of files

    Returns
    -------
    start_index : list
        indicate the index for the corresponding lead time (start)
    end_index : list
        indicate the index for the corresponding lead time (end)

    """
    compt = 0
    start_index = []
    end_index   = []
    
    for file in RDPS_files:
        # name of the fie
        name_rdps = nameprodfile['RDPS'][compt]
        
        # boolean list check if the required time window is the files
        lst = [x in file.variables["time"][:] for x in required_time]
        
        if(all(lst)):
            pass
        else:
            missinghours = list(np.where(~np.array(lst))[0])
            sys.exit(f"""RDPS file {name_rdps} must contain at least {len(required_time)} hours
                     from {required_time[0]} to {required_time[-1]} forecast time window - 
                     the following hour(s) {missinghours} is(are) missing""")

        compt       = compt + 1
        
        start_index.append(int(np.where(file.variables["time"][:]==required_time[1])[0]))
        end_index.append(int(np.where(file.variables["time"][:]==required_time[-1])[0]) + 1)
        
    return start_index, end_index
    
 
def create_netCDF(source_file, start_date, output_directory):
    """
    
    Creates a new NetCDF file based on a source file (RDRS or RDPS) and set up its dimensions
    Parameters
    ----------
    source_file : netcdf dataset
        can be either RDPS or RDRS.
    start_date : str
        string datasets.
    output_directory : str
        path of the output directory.

    Returns
    -------
    ncid : nc id
        identifiant of the netcdf file.

    """
    # get the desired name for the output NetCDF file
    start_date_str   = "".join(( start_date.strftime("%Y%m%d") , "12.nc"))                               
    output_file_name = "".join(("RDPS_CaPA_", start_date_str))
    
    # first if the file already exist delete it -- otherwise problem for creating netcdf with the same name
    outfile_nc       = join(output_directory, output_file_name)
    if exists(outfile_nc):
        os.remove(outfile_nc)

    #Check if output directory exist
    chech_f = os.path.isdir(out_directory)

    # Create a folder if doesn't exist
    if not chech_f:
        os.makedirs(output_directory)
        print("folder created : ", out_directory)

    else:
        print(out_directory, "folder already exists.")

        
    # create new netCDF file for output
    ncid             = nc4.Dataset(outfile_nc, "w", format="NETCDF4")
    
    # create dimensions for netCDF file
    for name, dimension in source_file.dimensions.items():
        ncid.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
        
    # set file-level attributes
    ncid.setncatts(source_file.__dict__)
    ncid.product = "RDPS_CaPA"
    print (str(source_file.__dict__['product']))
    
    # needs to be adopted for     RDRS_v2.1 of the reanalysis (now works only for RDRS_v2
    ncid.Remarks = "Variable names are following the convention <Product>_<Type:A=Analysis,P=Prediction>_<ECCC name>_<Level/Tile/Category>. Variables with level \'10000\' are at surface level. The height [m] of variables with level \'0XXXX\' needs to be inferrred using the corresponding fields of geopotential height (GZ_0XXXX-GZ_10000). The variables UUC, VVC, UVC, and WDC are not modelled but inferred from UU and VV for convenience of the users. Precipitation (PR) is reported as 6-hr accumulations for CaPA_fine and CaPA_coarse. Precipitation (PR) are accumulations since beginning of the forecast for GEPS, GDPS, REPS, RDPS, HRDPS, and CaLDAS. The re-analysis product RDRS (_v1, _v2, v2.1) contains two variables for precipitation: \'P_PR_SFC\' is the model precipitation (trial field used by CaPA) and \'A_PR_SFC\' is precipitations adjusted with CaPA 24h precipitation. Please be aware that the baseflow \'O1\' of the current version of WCPS is not reliable during the spring melt period. RDPS_CaPA is not a product available directly on CaSPar. It is created by combining RDPS and CaPA outputs."
    
    return ncid

def netCDF_variable_assignment(ncid, source_file, var_list, RDPS_files, RDRS_file, interpolate_bool, start_index, end_index):
    """
    Creates the variables in the output NetCDF file
    Parameters
    ----------
    ncid : nc id
        identifiant of the netcdf file.
    source_file : netcdf dataset
        can be either RDPS or RDRS.
    var_list : list
        list of variables to be added in the file output.
    RDPS_files : list
        list of RDPS netcdf datasets.
    RDRS_file : netcdf dataset
        RDRS template used if interpolation is needed.
    interpolate_bool : bool
        if true all variables (from RDPS) will be interpolated to the RDRS grid.
    start_index : list
        list of the start index (to take the right time step - start index).
    end_index : list
        list of the start index (to take the right time step - end index).

    Returns
    -------
    None.

    """
    # create time variable - metadata has to be from the RDRS file
    # the time dimension is unlimited
    times = ncid.createVariable('time', 'i4', ('time',))
    times.setncatts(RDRS_file['time'].__dict__)
    times[:] = np.arange(1, 25, 1)       
    
    # create variable holding location of the rotated pole
    pole = ncid.createVariable('rotated_pole', 'f4')
    pole.setncatts(source_file['rotated_pole'].__dict__)
    
    # create spatial variables and attributes
    for latlon in ('lat','lon'):
        rlatlon          = 'r'+latlon
        ncid.createVariable(rlatlon,'f4',(rlatlon,))
        
        ncid[rlatlon].setncatts(source_file[rlatlon].__dict__)
        ncid[rlatlon][:] = source_file[rlatlon][:]
        
        ncid.createVariable(latlon,'f4',('rlat','rlon'))
        ncid[latlon].setncatts(source_file[latlon].__dict__)
        ncid[latlon][:] = source_file[latlon][:]

    # create all other variables (except precipitation analysis)
    # loop through each variable name from input RDPS_files
    for var in RDPS_files[0].variables:
        if var in var_list:       
            variable_data(ncid, var, RDPS_files, RDRS_file, interpolate_bool, start_index, end_index)

            
def lon_lat_to_cartesian(lon, lat):
    """
    
    Parameters
    ----------
    lon : list
        longitude values.
    lat : list
        latitude values.

    Returns
    -------
    x : list
        coordinates in euclidian dist.
    y : list
        coordinates in euclidian dist..
    z : list
        coordinates in euclidian dist..

    """
    # WGS 84 reference coordinate system parameters
    A   = 6378.137 # major axis [km]   
    E2  = 6.69437999014e-3 # eccentricity squared 
    
    lon_rad = np.radians(lon)
    lat_rad = np.radians(lat)
    
    # convert to cartesian coordinates
    r_n = A / (np.sqrt(1 - E2 * (np.sin(lat_rad) ** 2)))
    x = r_n * np.cos(lat_rad) * np.cos(lon_rad)
    y = r_n * np.cos(lat_rad) * np.sin(lon_rad)
    z = r_n * (1 - E2) * np.sin(lat_rad)
    return x,y,z        

def interpolate_grid1_to_grid2(lat_grid1, lon_grid1, lat_grid2, lon_grid2):
    """
    Interpolation using nearest neighbor
    Parameters
    ----------
    lat_grid1 : list or np.array
        latitude of the grid 1 (initial grid).
    lon_grid1 : list or np.array
        longitude of the grid 1 (initial grid).
    lat_grid2 : list or np.array
        latitude of the grid 2 (final grid).
    lon_grid2 : list or np.array
        longitude of the grid 2 (final grid)..

    Returns
    -------
    d : list
        distance between points from grid1 to grid2.
    inds : list
        list of index used to identify the right index.

    """
    
    # Make sure longitudes are between 0 and 360
    lon_grid1     = lon_grid1 % 360
    lon_grid2     = lon_grid2 % 360
    
    # transform into cartesian coordinates
    x1, y1, z1  = lon_lat_to_cartesian(lon_grid1.flatten(), lat_grid1.flatten())
    x2, y2, z2  = lon_lat_to_cartesian(lon_grid2.flatten(), lat_grid2.flatten())
    
    # build the kdtree
    tree        = KDTree(np.column_stack((x1, y1, z1)))
    
    # build the query based on the nearest neighbor
    d, inds     = tree.query(np.column_stack((x2, y2, z2)), k = 1)
    
    return d, inds
    
def variable_data(ncid, key_var, RDPS_files, RDRS_file, interpolate_bool, start_index, end_index):

    """
    Assigns data to the variables of the output netcdf file
    Parameters
    ----------
    ncid : nc id
        identifiant of the netcdf file.
    key_var : str
        name of the variable to be added in the file output.
    RDPS_files : list
        list of RDPS netcdf datasets.
    RDRS_file : netcdf dataset
        RDRS template used if interpolation is needed.
    interpolate_bool : bool
        if true all variables (from RDPS) will be interpolated to the RDRS grid.
    start_index : list
        list of the start index (to take the right time step - start index).
    end_index : list
        list of the start index (to take the right time step - end index).

    Returns
    -------
    None.
    """
    
    # create variable for current key (dummy variable)
    ncid.createVariable("RDPS_CaPA_"+key_var[5:], 'f4', ('time', 'rlat', 'rlon'))
    
    # set attributes of variable
    ncid["RDPS_CaPA_"+key_var[5:]].setncatts(RDPS_files[0][key_var].__dict__)
    
    # create base array of data from first RDPS file
    data = np.array(RDPS_files[0].variables[key_var][:][start_index[0] : end_index[0]])
    
    # get full data
    for n in range(len(RDPS_files)):
        if n != 0:
            data = np.append(data, np.array(RDPS_files[n].variables[key_var][:][start_index[n] : end_index[n]]), axis = 0)
    
    # change grid system to the RDRS one
    if interpolate_bool:
        
        logging.info(f"Doing the interpolation from RDPS to RDRS grid for variable: {key_var}")
        # initial grid
        rdps_lon = RDPS_files[0].variables['lon'][:]
        rdps_lat = RDPS_files[0].variables['lat'][:]
        # destination grid
        rdrs_lon = RDRS_file.variables['lon'][:]
        rdrs_lat = RDRS_file.variables['lat'][:]
        
        # identify the right index
        dist, idx = interpolate_grid1_to_grid2(rdps_lon, rdps_lat, rdrs_lon, rdrs_lat)
    
        # do the interpolation
        data_grid_rdrs = []
        
        for n in range(len(data)):
            
            data_int = data[n].flatten()[idx]
            
            data_int2 = np.where(dist>15, np.nan, data_int)

            data_0    = data_int2.reshape(rdrs_lon.shape)
        
            data_grid_rdrs.append(data_0)
    
        
    ncid["RDPS_CaPA_"+key_var[5:]][:] = data_grid_rdrs        
    
       
def do_A_PR_SFC(ncid, RDPS_files, CaPA_files, RDRS_file, interpolate_bool):
    """
    Create the "Analysis: Quantity of precipitation" variable in the output netCDF file
    Parameters
    ----------
    ncid : nc id
        identifiant of the netcdf file.
    RDPS_files : list
        list of RDPS netcdf datasets.
    CaPA_files : list
        list of CaPA netcdf datasets.
    interpolate_bool : bool
        if true all variables (from RDPS) will be interpolated to the RDRS grid.

    Returns
    -------
    None.
    """
    
    # number of hour
    number_of_hours = 6
    # number of files in list
    number_of_files = len(CaPA_files)
    
    # create a list for the RDPS precipitation values
    rdps_data       = []
    for file in RDPS_files:
        rdps_data.append(file.variables['RDPS_P_PR_SFC'][:])
    
    # match dimensions of the CaPA and RDPS data by using interpolation
    # set-up variables for interpolation
    capa_lon    = CaPA_files[0].variables['lon'][:]
    capa_lat    = CaPA_files[0].variables['lat'][:]
    
    rdps_lon    = RDPS_files[0].variables['lon'][:]
    rdps_lat    = RDPS_files[0].variables['lat'][:]
    
    # identify the right index
    dist, idx   = interpolate_grid1_to_grid2(capa_lon, capa_lat, rdps_lon, rdps_lat)
    
    # create a list of the CaPA variables
    logging.info('BEGIN the interpolation of CaPA precipitation to RDPS grid')
    capa_data   = []
    for n  in range(number_of_files):
        capa_pr          = CaPA_files[n].variables['CaPA_coarse_A_PR_SFC'][:][0]

        capa_to_rdpsgrid = capa_pr.flatten()[idx].reshape(rdps_lon.shape)

        capa_data.append(capa_to_rdpsgrid)
    logging.info('END the interpolation of CaPA precipitation to RDPS grid')
    
    # create empty arrays for the new precipitation data
    precip_data_hourly = []
    for item in rdps_data:
        precip_data_hourly.append(np.empty_like(item)[0:number_of_hours,:,:]) # ATTENTION PAS DE 0 À 6 
    
    # get hourly precipitation forecasts
    for n in range(number_of_files):
        for i in range(number_of_hours):
            precip_data_hourly[n][i] = rdps_data[n][i + 1] - rdps_data[n][i]
    
    # get RDPS forecasted total precipitation for the same 6 hour window as capa
    rdps_precip_daily_sums = []
    for item in precip_data_hourly:
        rdps_precip_daily_sums.append(np.sum(item, axis = 0))

    # re-scale hourly RDPS forecasted total precipitation using a factor based on
    # the 6-h accumulation of both capa and rdps
    for n in range(number_of_files):
        
        scale       = capa_data[n]/rdps_precip_daily_sums[n]
        capa_hr_cst = capa_data[n]/number_of_hours
        
        for t in range(number_of_hours):
            precip_data_hourly[n][t] = np.where(rdps_precip_daily_sums[n]>1e-6, scale*precip_data_hourly[n][t], capa_hr_cst)    
            
                       
    precip_data = np.array(precip_data_hourly)
    precip_data = np.concatenate(precip_data, axis = 0)
    
    # interpolation of the hourly precipitation onto the RDRS grid
    if interpolate_bool:
        
        # points for interpolation
        rdps_lon = RDPS_files[0].variables['lon'][:]
        rdps_lat = RDPS_files[0].variables['lat'][:]
        rdrs_lon = RDRS_file.variables['lon'][:]
        rdrs_lat = RDRS_file.variables['lat'][:]
        
        # identify the right index
        dist, idx   = interpolate_grid1_to_grid2(rdps_lon, rdps_lat, rdrs_lon, rdrs_lat)
               
        # do the interpolation
        pr_hr_grid_rdrs = []
        
        for n in range(len(precip_data)):
            
            precip_data_int     = precip_data[n].flatten()[idx]
            
            precip_data_int2    = np.where(dist>15, np.nan, precip_data_int)

            pr_hr_grid_rdrs_0   = precip_data_int2.reshape(rdrs_lon.shape)
            
            # pr_hr_grid_rdrs_0 = precip_data[n].flatten()[idx].reshape(rdrs_lon.shape)
        
            pr_hr_grid_rdrs.append(pr_hr_grid_rdrs_0)
        
        precip_data = pr_hr_grid_rdrs
    
    # create precipitation variable
    ncid.createVariable('RDPS_CaPA_A_PR_SFC', 'f4', ('time', 'rlat', 'rlon'))
    
    name_product=RDRS_file.__dict__['product']
    name_product_PR=str(name_product)+'_A_PR0_SFC'
       
    #ncid['RDPS_CaPA_A_PR_SFC'].setncatts(RDRS_file['RDRS_v2_A_PR0_SFC'].__dict__)
    ncid['RDPS_CaPA_A_PR_SFC'].setncatts(RDRS_file[name_product_PR].__dict__)
    ncid['RDPS_CaPA_A_PR_SFC'][:] = precip_data
    
    
    return(precip_data)

def do_P_TD(ncid, RDPS_files, CaPA_files, RDRS_file, interpolate_bool, start_index, end_index):

    """
    Create the "Forecast: Dew point temperature 09950" variable in the output netCDF file
    Parameters
    ----------
    ncid : nc id
        identifiant of the netcdf file.
    RDPS_files : list
        list of RDPS netcdf datasets.
    CaPA_files : list
        list of CaPA netcdf datasets.
    interpolate_bool : bool
        if true all variables (from RDPS) will be interpolated to the RDRS grid.
    start_index : list
        list of the start index (to take the right time step - start index).
    end_index : list
        list of the start index (to take the right time step - end index).

    Returns
    -------
    None.
    """
    # create base array of data from first RDPS file for the temperature
    temp_data = np.array(RDPS_files[0].variables["RDPS_P_TT_09950"][:][start_index[0] : end_index[0]])
    
    # create base array of data from first RDPS file for the relative humidity
    RH_data   = np.array(RDPS_files[0].variables["RDPS_P_HR_09950"][:][start_index[0] : end_index[0]])
    
    # get full data
    for n in range(len(RDPS_files)):
        if n != 0:
            temp_data = np.append(temp_data, np.array(RDPS_files[n].variables["RDPS_P_TT_09950"][:][start_index[n] : end_index[n]]), axis = 0)
            RH_data = np.append(RH_data, np.array(RDPS_files[n].variables["RDPS_P_HR_09950"][:][start_index[n] : end_index[n]]), axis = 0)
            
    # calculate saturation vapor pressure in mb
    es = 6.112 * (math.e ** ((17.67 * temp_data) / (temp_data + 243.5)))
    
    # calculatre vapor pressure in mb
    e = es * RH_data
    
    # calculate Dew Point Temperature in deg C
    TD_data = np.log(e / 6.112) * 243.5 / (17.67 - np.log(e / 6.112))
    
    # interpolation
    if interpolate_bool:
        # points for interpolation
        rdps_lon = RDPS_files[0].variables['lon'][:]
        rdps_lat = RDPS_files[0].variables['lat'][:]
        rdrs_lon = RDRS_file.variables['lon'][:]
        rdrs_lat = RDRS_file.variables['lat'][:]
        
        # identify the right index
        dist, idx   = interpolate_grid1_to_grid2(rdps_lon, rdps_lat, rdrs_lon, rdrs_lat)
        
        
        # do the interpolation
        TD_data_rdrs_grid = []
        for n in range(len(TD_data)):
            TD_data_int = TD_data[n].flatten()[idx]
            TD_data_int2 = np.where(dist>15, np.nan, TD_data_int)
            TD_data_rdrs_grid_0 = TD_data_int2.reshape(rdrs_lon.shape)
            
            #TD_data_rdrs_grid_0 = TD_data[n].flatten()[idx].reshape(rdrs_lon.shape)
        
            TD_data_rdrs_grid.append(TD_data_rdrs_grid_0)

        # for n in range(len(data)):
            
        #     data_int = data[n].flatten()[idx]
            
        #     data_int2 = np.where(dist>15, np.nan, data_int)

        #     data_0    = data_int2.reshape(rdrs_lon.shape)
        
        #     data_grid_rdrs.append(data_0)


        TD_data = TD_data_rdrs_grid
        
    ncid.createVariable('RDPS_CaPA_P_TD_09950', 'f4', ('time', 'rlat', 'rlon'))
        
    name_product=RDRS_file.__dict__['product']
    print (name_product)
    name_product_TD=str(name_product)+'_P_TD_09944'
    
    #ncid['RDPS_CaPA_P_TD_09950'].setncatts(RDRS_file['RDRS_v2_P_TD_09944'].__dict__)
    ncid['RDPS_CaPA_P_TD_09950'].setncatts(RDRS_file[name_product_TD].__dict__)
    ncid['RDPS_CaPA_P_TD_09950'][:] = TD_data


# def main():
#%%  READ THE CONFIGURATION
config_1 = ConfigParser()
config_1.read("Configuration.ini")

# interpolate_bool is True if using RDRS grid definition, and False if using RDPS grid definition
interpolate_bool  = config_1.getboolean("Settings",  "interpolate_bool")
if interpolate_bool: 
    logging.info("You will be using RDRS grid definition for your outputs \n For changes look at the configuration file")
else:
    logging.info("You will be using RDPS grid definition for your outputs \n For changes look at the configuration file")


# setting to include or not include the precipitation analysis variable
analysis_bool       = config_1.getboolean("Settings", "precipitation_analysis_bool")
if analysis_bool: 
    logging.info("The precipitation analysis will be included in the output files ")
else:
    logging.info("WARNING: The precipitation analysis will NOT be included in the output files ")
        
# setting to include extra dew point temperature variable
dew_bool            = config_1.getboolean("Settings","dewpoint_09950_bool")
if dew_bool: 
    logging.info("The dew point temperature variable at that level will be included in the output files")
else:
    logging.info("WARNING: he dew point temperature variable at that level will NOT be included in the output files ")

# paths to the input directories and files
RDPS_files_directory = config_1["Input-files"]["RDPS_files_directory"]
CaPA_files_directory = config_1["Input-files"]["CaPA_files_directory"]
RDRS_file_path       = config_1["Input-files"]["RDRS_file_path"]

# list of variables to include in the output NetCDF file
var_list             = config_1["Variables"]["variables"].split(', ')
     
# start and end values
start_year                = config_1["Time"]["start_year"]
start_month               = config_1["Time"]["start_month"]
start_day                 = config_1["Time"]["start_day"]
end_year                  = config_1["Time"]["end_year"]
end_month                 = config_1["Time"]["end_month"]
end_day                   = config_1["Time"]["end_day"]
required_rdps_time_window = [int(x)  for x in config_1["Time"]["required_rdps_time_window"].split(',')]

# get the path to the output directory
output_directory = config_1["Output"]["output_directory"]
#%% 
# start and end date
start_date          = date(int(start_year),int(start_month),int(start_day))
end_date            = date(int(end_year),int(end_month),int(end_day))
 

while start_date <= end_date:
    
    # read the input netCDF files
    RDPS_files, CaPA_files, RDRS_file_tmp, namefiledict = read_data(RDPS_files_directory, CaPA_files_directory, RDRS_file_path, start_date)
    
    # test to see if all variables are present in the RDPS file
    check_vars(RDPS_files, var_list, namefiledict)
    
    # test to see if all RDPS files and CaPA files are on the same grid
    check_grid(RDPS_files, namefiledict, 'RDPS')
    check_grid(CaPA_files, namefiledict, 'CaPA') 
 
    # check to see if input RDPS files contain the required 7 hours and RDPS files contain the same forecast horizons
    # also figure out what array slices we need if file contains more than the required 7 hours

    start_index, end_index = check_time(RDPS_files, required_rdps_time_window, namefiledict)

    # set the file to get the grid definition from
    if interpolate_bool :
        source_file = RDRS_file_tmp
    else:
        source_file = RDPS_files[0]

    # create output netCDF file
    ncid = create_netCDF(source_file, start_date, output_directory)

    # create variables in output netCDF file
    netCDF_variable_assignment(ncid, source_file, var_list, RDPS_files, RDRS_file_tmp, interpolate_bool, start_index, end_index)

    if analysis_bool:
        # create "Analysis: Quantity of precipitation" variable
        do_A_PR_SFC(ncid, RDPS_files, CaPA_files, RDRS_file_tmp, interpolate_bool)
        
    if dew_bool:
        # create "Forecast: Dew point temperature_09950 variable"
        do_P_TD(ncid, RDPS_files, CaPA_files, RDRS_file_tmp, interpolate_bool, start_index, end_index)

    # close the netCDF file
    ncid.close()
    
    # increment the start date
    start_date = start_date + timedelta(days = 1)



# start = time.time()
# # run the script
# if __name__ == '__main__':
#     main()  
# end = time.time()
# print(end - start)
