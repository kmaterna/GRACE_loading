#!/usr/env/bin python
# Functions that compute PREM loading displacements. 
# In the same format as the C program I use for analytical solutions of the half-space. 
# Mode 1: Single station, complex time-dependent load. Results in time-dependent output.
# Mode 2: Given load, MULTIPLE stations. Results in an image or profile of deformation. 

import numpy as np
import haversine
import math
import helper_functions


def compute_prem_load(station_array, load_array, disk_radius, greensfunc, max_distance):
    # For a given station and a given set of loads,
    # Compute the loading on a PREM earth structure
    # Max distance is in degrees
    total_disp_x = [];
    total_disp_y = [];
    total_disp_v = [];
    for i in range(len(station_array.lon)):  # for each station...

        # Take a set of rectangular loads and tile them into disk loads.
        CircleLoads = get_disks(load_array, disk_radius, max_distance, station_array.lon[i], station_array.lat[i]);
        print("Number of circles: %d" % len(CircleLoads.lon));

        disk_sum_v = 0;
        disk_sum_x = 0;
        disk_sum_y = 0;
        for j in range(np.size(CircleLoads.lon)):
            dist_from_center = haversine.distance([CircleLoads.lat[j], CircleLoads.lon[j]],
                                                  [station_array.lat[i], station_array.lon[i]]);  # reported in km
            azimuth = haversine.calculate_initial_compass_bearing((station_array.lat[i], station_array.lon[i]),
                                                                  (CircleLoads.lat[j], CircleLoads.lon[j]));
            myindex = find_closest_greens(greensfunc.km,
                                          dist_from_center);  # find the nearest answer in the Green's functions
            vdisp = greensfunc.vload[myindex];  # in m
            hdisp = greensfunc.hload[myindex];  # in m
            [xdisp, ydisp] = az_decompose(hdisp, (station_array.lat[i], station_array.lon[i]), (
                CircleLoads.lat[j], CircleLoads.lon[j]));  # Decompose radial into x and y components

            disk_sum_v = disk_sum_v + vdisp * CircleLoads.pressure[
                j] * 1000 / 9.81;  # THIS SHOULD BE DIVIDED BY 9.81! the PREM code already has it. Result in mm
            disk_sum_x = disk_sum_x + xdisp * CircleLoads.pressure[
                j] * 1000 / 9.81;  # THIS SHOULD BE DIVIDED BY 9.81! the PREM code already has it. Result in mm
            disk_sum_y = disk_sum_y + ydisp * CircleLoads.pressure[
                j] * 1000 / 9.81;  # THIS SHOULD BE DIVIDED BY 9.81! the PREM code already has it. Result in mm

        total_disp_x.append(disk_sum_x);
        total_disp_y.append(disk_sum_y);
        total_disp_v.append(disk_sum_v);
    return total_disp_x, total_disp_y, total_disp_v;  # should exactly match length of input station_array fields.


def get_disks(load_array, disk_radius, threshold_distance, sta_lon, sta_lat):
    # threshold_distance is in degrees
    circle_lons = [];
    circle_lats = [];
    circle_radius = [];
    circle_pressure = [];
    threshold_distance = threshold_distance * (np.pi / 180) * 6370;  # converting to km  (6370km = 1 radian).
    circle_area = np.pi * disk_radius * disk_radius;  # Area = pi * r^2
    square_side_length = np.sqrt(circle_area);  # in same units as disk_radius

    for i in range(len(load_array.lon)):  # for each rectangular load cell
        if abs(load_array.lat[i]) > 85:  # ignoring the arctic and antarctic
            continue;
        dist_from_center = haversine.distance([load_array.lat[i], load_array.lon[i]],
                                              [sta_lat, sta_lon]);  # reported in km
        if dist_from_center < threshold_distance:
            # if the load cell is less than a certain distance from the station, discretize it.

            center = [load_array.lon[i], load_array.lat[i]];
            load_akm = load_array.east_width[i] * 111.0 * np.cos(center[1] * np.pi / 180);
            load_bkm = load_array.north_width[i] * 111.0;

            # The placement of each disk-shaped tile within the load
            n_tiles_x = int(np.floor(load_akm * 2 / square_side_length));
            n_tiles_y = int(np.floor(
                load_bkm * 2 / square_side_length));  # number of tiles that fit in x and y dimensions of rectangle
            x_array = [];
            y_array = [];
            x_nodes = np.linspace(center[0] - load_array.east_width[i], center[0] + load_array.east_width[i],
                                  n_tiles_x + 1, endpoint=True);
            y_nodes = np.linspace(center[1] - load_array.north_width[i], center[1] + load_array.north_width[i],
                                  n_tiles_y + 1, endpoint=True);

            for j in range(len(x_nodes) - 1):
                x_array.append((x_nodes[j] + x_nodes[j + 1]) * 0.5);  # average of two nodes
            for j in range(len(y_nodes) - 1):
                y_array.append((y_nodes[j] + y_nodes[j + 1]) * 0.5);  # average of two nodes

            [centerx, centery] = np.meshgrid(x_array, y_array);
            centerx = np.reshape(centerx, (np.size(centerx), 1))
            centery = np.reshape(centery, (np.size(centery), 1))

            # HERE WE RE-PACKAGE INTO CIRCULAR LOADS
            for j in range(np.size(centerx)):
                circle_lons.append(centerx[j]);
                circle_lats.append(centery[j]);
                circle_radius.append(disk_radius);
                circle_pressure.append(load_array.pressure[i]);

    CircleLoads = helper_functions.Load_Array(lon=circle_lons, lat=circle_lats, east_width=circle_radius,
                                              north_width=[], pressure=circle_pressure);
    return CircleLoads;


def find_closest_greens(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
        return idx - 1;
    else:
        return idx;


def az_decompose(hdisp, station_pos, disk_pos):
    # both station_pos and disk_pos are tuples
    compass_bearing = haversine.calculate_initial_compass_bearing(station_pos, disk_pos);
    azimuth = 90 - compass_bearing;
    xdisp = hdisp * -np.cos(azimuth * np.pi / 180);
    ydisp = hdisp * -np.sin(azimuth * np.pi / 180);
    return [xdisp, ydisp];
