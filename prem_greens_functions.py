#!/usr/env/bin python
# Functions that compute PREM loading displacements. 
# In the same format as the C program I use for analytical solutions of the half-space. 
# Mode 1: Single station, complex time-dependent load. Results in time-dependent output.
# Mode 2: Given load, MULTIPLE stations. Results in an image or profile of deformation. 

import numpy as np 
import haversine
import math


def compute_prem_load(station_file, rect_loadfile, disk_loadfile, disk_radius, PREM_green, disk_out_file, threshhold_distance, mode):
	# For a given station and a given set of disk loads,
	# Compute the loading on a PREM earth structure

	# CONFIG AND INPUTS
	[sta_lon, sta_lat, dep] = np.loadtxt(station_file,unpack=True);
	[deg, km, vload, hload] = np.loadtxt(PREM_green,unpack=True);

	if mode==1:
		ofile=open(disk_out_file,'a');  # compute a time series for a given point
	if mode==2:
		ofile=open(disk_out_file,'w');  # compute an image at many points 

	# COMPUTE
	for i in range(np.size(sta_lon)):  # for each station... 

		# Logic for looping in case there's only one station
		if np.size(sta_lon)==1:
			sta_lon=[sta_lon]; sta_lat=[sta_lat]; dep=[dep];		
		
		# Take a set of rectangular loads and tile them into disk loads. 
		get_tiles_write_disk_loadfile(rect_loadfile, disk_loadfile, disk_radius, threshhold_distance, sta_lon[i], sta_lat[i]);  # set up tiles

		# Logic for looping in case there's only one disk
		[disk_lon, disk_lat, radius, p] = np.loadtxt(disk_loadfile,unpack=True);
		if np.size(disk_lon)==1:
			disk_lon=[disk_lon]; disk_lat=[disk_lat]; radius=[radius]; p=[p];

		disk_sum_v=0; disk_sum_x=0; disk_sum_y=0;
		
		for j in range(np.size(disk_lon)):

			dist_from_center=haversine.distance([disk_lat[j], disk_lon[j]],[float(sta_lat[i]), float(sta_lon[i])]);  # reported in km
			azimuth=haversine.calculate_initial_compass_bearing((float(sta_lat[i]), float(sta_lon[i])), (disk_lat[j], disk_lon[j]))
			
			myindex = find_closest_greens(km, dist_from_center);   # find the nearest answer in the Green's functions
			vdisp=1000 * vload[myindex];  # in mm
			hdisp=1000 * hload[myindex];  # in mm
			[xdisp, ydisp] = az_decompose(hdisp,(float(sta_lat[i]), float(sta_lon[i])), (disk_lat[j], disk_lon[j])); # Decompose radial into x and y components
			
			disk_sum_v=disk_sum_v+vdisp*p[j]/9.81; # THIS SHOULD BE DIVIDED BY 9.81! the PREM code already has it. 
			disk_sum_x=disk_sum_x+xdisp*p[j]/9.81; # THIS SHOULD BE DIVIDED BY 9.81! the PREM code already has it. 
			disk_sum_y=disk_sum_y+ydisp*p[j]/9.81;# THIS SHOULD BE DIVIDED BY 9.81! the PREM code already has it. 

		ofile.write("%f %f %f %f %f %f\n" % (sta_lon[i], sta_lat[i], dep[i], disk_sum_x, disk_sum_y, disk_sum_v) );
	ofile.close();
	return;



def get_tiles_write_disk_loadfile(rect_loadfile, disk_loadfile, disk_radius, threshhold_distance, station_lon, station_lat):
	# Output Format = [lon lat radius pressure]
	# This takes a series of rectangular loads and computes the equivalent disk loads. 
	# Tile an area with a bunch of circles with area that matches square tiles. 
	
	ofile=open(disk_loadfile,'w');

	[load_lon, load_lat, a_deg, b_deg, p] = np.loadtxt(rect_loadfile,unpack=True);  # one line for each rectangle you want to tile/write
	
	if np.size(load_lon)==1:
		load_lon=[load_lon]; load_lat=[load_lat]; a_deg=[a_deg]; b_deg=[b_deg]; p=[p]; 

	for i in range(len(load_lon)): # for each rectangular load cell
		dist_from_center=1000*haversine.distance([load_lat[i], load_lon[i]],[float(station_lat), float(station_lon)]);  # reported in km, so we multiply to meters
		if dist_from_center < threshhold_distance:  # if the load cell is less than a certain distance from the station, discretize it. 
	
			circle_area = np.pi*disk_radius*disk_radius;  # Area = pi * r^2  ( in m)
			square_side_length = np.sqrt(circle_area);

			center=[load_lon[i], load_lat[i]];
			load_akm = a_deg[i]*111.0*np.cos(center[1]*np.pi/180);
			load_bkm = b_deg[i]*111.0;

			# The placement of each disk-shaped tile within the load
			n_tiles_x = np.floor(load_akm*2 / square_side_length);
			n_tiles_y = np.floor(load_bkm*2 / square_side_length);  # the number of tiles that fit in x and y dimensions of the rectangle	
			x_array=[]; y_array=[];
			x_nodes = np.linspace(center[0]-a_deg[i], center[0]+a_deg[i], n_tiles_x+1, endpoint=True);
			y_nodes = np.linspace(center[1]-b_deg[i], center[1]+b_deg[i], n_tiles_y+1, endpoint=True);

			for j in range(len(x_nodes)-1):
				x_array.append((x_nodes[j]+x_nodes[j+1])*0.5);  # average of two nodes
			for j in range(len(y_nodes)-1):
				y_array.append((y_nodes[j]+y_nodes[j+1])*0.5);  # average of two nodes

			[centerx,centery]=np.meshgrid(x_array, y_array);
			centerx=np.reshape(centerx,(np.size(centerx),1))
			centery=np.reshape(centery,(np.size(centery),1))
		
			for j in range(np.size(centerx)):
				ofile.write("%f %f %f %f\n" % (centerx[j], centery[j], disk_radius, p[i]) );   

	ofile.close();
	return;


def find_closest_greens(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1;
    else:
        return idx;

def az_decompose(hdisp,station_pos, disk_pos):
	# both station_pos and disk_pos are tuples
	compass_bearing = haversine.calculate_initial_compass_bearing(station_pos, disk_pos);
	azimuth = 90-compass_bearing;
	xdisp=hdisp*-np.cos(azimuth*np.pi/180);
	ydisp=hdisp*-np.sin(azimuth*np.pi/180);
	return [xdisp, ydisp];


def km2deg_xdistance(km_distance, lat):
	# How many degrees is the x-distance in km?
	deg=(km_distance/111.0)/np.cos(lat*np.pi/180);
	return deg;

def km2deg_ydistance(km_distance):
	# How many degrees is the y-distance in km?
	deg=(km_distance/111.0);
	return deg;
