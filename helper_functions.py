#!/usr/env/bin python

# September 18, 2017
# Helper functions for the elastic loading time series calculations
# These are common to all calculations, no matter what Green's functions we use. 


import datetime as dt
from scipy.io import netcdf
import collections
from subprocess import call
import glob
import numpy as np

Station_Array = collections.namedtuple('Station_Array', ['name', 'lon', 'lat', 'start_time', 'end_time']);
Load_Array = collections.namedtuple('Load_Array',
                                    ['lon', 'lat', 'east_width', 'north_width', 'pressure']);  # for rectangular loads
GRACE_Array = collections.namedtuple('GRACE_Array', ['grace_central_times',
                                                     'layerstrings', 'x_range', 'y_range', 'we', 'east_width',
                                                     'north_width', 'scalefactor']);
Greensfunc = collections.namedtuple('Greensfunc', ['deg', 'km', 'hload', 'vload']);


def read_GRACE_netfile(netfile, scalefile):
    print("Reading GRACE file %s " % netfile);
    f = netcdf.netcdf_file(netfile, 'r');
    x_range = f.variables['lon'];  # an array from 0-360 in half-degree increments
    x_range = x_range[:].copy();
    y_range = f.variables['lat'];  # an array from -89 to 89 in half-degree increments
    y_range = y_range[:].copy();
    z_range = f.variables['time'];  # an array of 182 numbers, units = days since 2002-01-01
    z_range = z_range[:].copy();
    official_start_date = dt.datetime.strptime("2002-01-01", "%Y-%m-%d");
    zdates = [];
    layerstrings = [];
    for i in range(len(z_range)):
        zdates.append(official_start_date + dt.timedelta(days=int(z_range[i])));
        layerstrings.append(dt.datetime.strftime(zdates[-1], "%Y-%m-%d"));
    # the actual middle date of each interval, matching with the excel file for grace months

    we = f.variables['lwe_thickness'];
    we = we[:].copy();  # shape 182 x 360 x 720.  360 for lat, 720 for lon, 182 for time.
    north_interval = 180 / np.shape(we)[1];
    east_interval = 360 / np.shape(we)[2];
    east_width = east_interval / 2 * np.ones(np.shape(we[0, :, :]));  # in degrees, half-width of the box
    north_width = north_interval / 2 * np.ones(np.shape(we[0, :, :]));  # in degrees, half-width of the box

    # Get the GRACE scale factor corrections.
    if "Mascons" in netfile:
        scalefactor = [];
    else:
        sf = netcdf.netcdf_file(scalefile, 'r');
        scalefactor = sf.variables['SCALE_FACTOR'];  # size = ; time-independent.
        scalefactor = scalefactor[:].copy();

    grace_data_structure = GRACE_Array(grace_central_times=zdates, layerstrings=layerstrings,
                                       x_range=x_range, y_range=y_range, we=we, east_width=east_width,
                                       north_width=north_width, scalefactor=scalefactor);
    return grace_data_structure;


def read_station_file(station_list):
    print("Reading station file %s " % station_list);
    ifile = open(station_list, 'r');
    station_names = [];
    station_lon = [];
    station_lat = [];
    start_time = [];
    end_time = [];
    grace_begin = dt.datetime.strptime("2002-01-01", "%Y-%m-%d");
    grace_end = dt.datetime.strptime("2017-06-01", "%Y-%m-%d");
    grace_default_start = dt.datetime.strptime("2012-01-01", "%Y-%m-%d");
    for line in ifile:
        temp = line.split();
        station_names.append(temp[0]);
        station_lon.append(float(temp[1]));
        station_lat.append(float(temp[2]));

        # Decide on the time range for the computation at each GPS station.
        if len(temp) > 4:
            gps_start_time = dt.datetime.strptime(temp[3], "%Y-%m-%d");
            gps_end_time = dt.datetime.strptime(temp[4], "%Y-%m-%d");
            gps_start_time = gps_start_time - dt.timedelta(days=400);  # giving some grace period
            gps_end_time = gps_end_time + dt.timedelta(days=400);

            # We can't compute GRACE before 2002 or after 2017, so ignore all data before that.
            if gps_start_time < grace_begin:
                gps_start_time = grace_begin;
            if gps_end_time > grace_end:
                gps_end_time = grace_end;

        else:  # If time constraints are not given, replace them with default time window.
            gps_start_time = grace_default_start;
            gps_end_time = grace_end;

        start_time.append(gps_start_time);
        end_time.append(gps_end_time);
        my_station_array = Station_Array(name=station_names, lon=station_lon,
                                         lat=station_lat, start_time=start_time, end_time=end_time);
    return my_station_array;


def read_spherical_greensfunc(gf_file):
    [deg, km, vload, hload] = np.loadtxt(gf_file, unpack=True);
    GF = Greensfunc(deg=deg, km=km, hload=hload, vload=vload);
    return GF;


def package_loads_from_GRACE(grace_data_structure, n, use_scale):
    """
    This reads the n'th slice from the netfile [X x Y x n], multiplies by the scale file,
    and packages the load into a basic object.
    If use_scale is set, then we multiply the GRACE timeslice by the scale factor grid.
    This is a global load.
    The widths are in degrees.
    """
    snapshot = grace_data_structure.we[n, :, :];  # Get a snapshot
    s1 = np.shape(snapshot);  # fill in NaN's as necessary over the water mask.
    mult_scalefactor = grace_data_structure.scalefactor;
    if use_scale == 0:  # If no scaling grid applied:
        mult_scalefactor = np.ones(s1);  # if we are not using scale factor, we just multiply everything by 1.

    load_lon = [];
    load_lat = [];
    load_width = [];
    load_height = [];
    load_pressure = [];
    for i in range(s1[0]):
        for j in range(s1[1]):
            if snapshot[i, j] < 5000:  # DOES THIS STILL WORK FOR MASCONS?
                pressure = snapshot[i, j] * mult_scalefactor[i, j] * (1.0 / 100) * 9.81 * 1000.0;
                # rho*g*h (h converted cm --> meters)
                load_lon.append(grace_data_structure.x_range[j]);
                load_lat.append(grace_data_structure.y_range[i]);
                load_width.append(grace_data_structure.east_width[i, j]);
                load_height.append(grace_data_structure.north_width[i, j]);
                load_pressure.append(pressure);
    load_obj = Load_Array(lon=load_lon, lat=load_lat, east_width=load_width, north_width=load_height,
                          pressure=load_pressure);
    return load_obj;


def write_global_load_file(outfile, Load_Array):
    # This function might be useful eventually, but I don't have reason to use it right now.
    return;


def delete_files_matching(match_string):
    clear_list = glob.glob(match_string);
    for item in clear_list:
        call(['rm', item], shell=False);
    return;
