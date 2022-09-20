# Instructions on running python script

## Objective of the script

This script is meant to combine RDPS and CaPA outputs to mimic RDRS.
RDPS and CaPA files are downloaded from https://caspar-data.ca/caspar. Details are given bellow.
To use the script, you will need to use the Configuration.ini file to provide inputs.
You will need to enter the paths to the directory containing RDPS files, the directory containing CaPA files (CaPA_coarse from CaSPAr), and the path to RDRS file. 

The followings indicate the different items to take into account before running the python script.

**1. Python setup**

The python script runs with 3.8 version of python and needs the following librairies:
- [numpy](https://numpy.org/)
- [scipy](https://scipy.org/)
- [netCDF4](https://pypi.org/project/netCDF4/)

**2. Download the needed data from CaSPAR**

As the objective is to build RDPS files that mimic RDRS files for different variables, you will need to download:
1. RDPS-Regional Deterministic Prediction System NetCDF files that contain forecast horizons of 6-12. Each forecast issue will be used for the forecast horizons 7-12 of that same time.
2. CaPA_coarse-Canadian Precipitation Analysis North America 10 km NetCDF files for all forecast issues and forecast horizons.

For example, if you need to build an RDPS that mimic RDRS for the date 2015010112.nc which contains hourly data from 2015-01-01 at 13UTC to 2015-01-02 at 12UTC to; you will need:
1. RDPS: 2015010106.nc, 2015010112.nc, 2015010118.nc and 2015010200.nc (all containing the 6 to 12 forecast horizons)
2. CaPA coarse: 2015010118.nc, 2015010200.nc, 2015010206.nc and 2015010212.nc
3. Optional - if you interpolate to RDRS grid (`interpolate_bool` is set to True); you will need RDRS file as template

**3. Adapt the configuration file**

You need to modify the following items in the configuration file:
1. Change the directory paths, which correspond to the places where you downloaded data from CaSPAR (see section above)
2. Change the `RDRS_file_path`, it is a file that is used as a template for the grid definition of your outputs. This file is used only if `interpolate_bool` is set to True. So if you want that your outputs follow the same grid as RDRS, you can choose an RDRS file that you already have (it can be the one of your choice for a given date) and set  `interpolate_bool=True`
3. In the `Time` section, you need to choose your first and last date of the period you want to create your outputs. The variable `required_rdps_time_window` does not need to be changed
4. The `variables` entry list of the variable that will be in your output files. For example, if  -- before year 2018 -- you had in your RDRS files the variables `RDPS_P_FB_SFC`, `RDPS_P_FI_SFC` and `RDPS_P_PR_SFC` then you will have `variables=RDPS_P_FB_SFC, RDPS_P_FI_SFC,RDPS_P_PR_SFC` if you want to keep them in your outputs for years after 2018. The only thing is to separate each variable by a comma. In the default configuration file all possible variables are provided. 
5. The settings for the interpolation are as follow:
    - If you want that your outputs to be interpolated on the RDRS grid, set `interpolate_bool=True`. Otherwise (`interpolate_bool=False`), your outputs will be on the same grid as your inputs, meaning the same as RDPS grid.
    - If you want to add the precipitation analysis (based on CaPA coarse) in your outputs, you need to set `precipitation_analysis_bool = True`. Otherwise (`precipitation_analysis_bool = False`), the precipitation analysis will not be added to your output file.
    - If you want to add the dew point temperature in your outputs file you will need to set (`dewpoint_09950_bool = True`). This variable is not available in RDPS but can be computed based on temperature and relative humidity from RDPS. If  `dewpoint_09950_bool = True`, be careful to have in your downloaded RDPS file the following variables:`RDPS_P_TT_09950` and `RDPS_P_HR_09950`. If you do not want the dew point temperature, set `dewpoint_09950_bool = False`.
6. Finally, the `output_directory` variable must be changed for the path were you want to put your output files.


**3. Example**

The folder `example` contains examples of RDRS, RDPS and CaPA files that can be used to generate RDPS that mimic RDRS for the 2015-01-01 date. The output will be generated in the folder `example/dataout`. 



