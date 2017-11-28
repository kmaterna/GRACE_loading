#!/usr/env/bin python

# EVALUATE GRACE LOADS ON GPS STATIONS, PREM STRUCTURE
# INPUT/network.txt contains list of: GPS_name, lon, lat, start_time, end_time
# For each timestamp after starting time, 
# We compute GRACE load at the GPS station from PREM earth structure
# timestamp   --> output
# Deformation --> output 
# Written by Kathryn Materna, 2017


import numpy as np 
from subprocess import call 
import glob
import datetime as dt 
import os
import matplotlib.pyplot as plt 
import prem_greens_functions
import helper_functions


def prem_earth_grace_timeseries(network, dataset, use_scale):
	[netfile, scalefile, station_file, output_dir, file_prefix, disk_radius, PREM_green, threshhold_distance] = configure(network, dataset, use_scale);
	[my_station_array, grace_data_structure]   = get_inputs(station_file, netfile, scalefile);
	compute(output_dir, file_prefix, use_scale, disk_radius, PREM_green, threshhold_distance, my_station_array, grace_data_structure);
	outputs(output_dir, file_prefix, my_station_array);
	return;



# ---------------------- CONFIGURE ------------------------ # 
def configure(network, dataset, use_scale):

	# GET STARTED WITH SET-UP AND CLEAN-UP
	input_file="INPUT/"+network+".txt"
	[netfile, scalefile, output_dir, file_prefix] = helper_functions.define_io_options(network, dataset, use_scale);
	helper_functions.delete_files_matching(output_dir+file_prefix+'*_PREM_model_ts.txt');  # get started with a clean directory. 

	# Other parameters that don't often change. 
	threshhold_distance = 2000*1000;  # in m. Further than 2000km away, we assume the elastic loads don't matter. 	
	disk_radius = 6.206;  # options are 0.564km (1km^2) or 31km (3019km^2), because that's what we have green's functions for. 
	PREM_green ="Wahr_6.206"  # Make sure these match!!!

	return [ netfile, scalefile, input_file, output_dir, file_prefix, disk_radius, PREM_green, threshhold_distance];



# ---------------------- INPUTS ------------------------ # 
def get_inputs(station_list, netfile, scalefile):
	# Bring us arrays of lat, lon, and time ranges that we need to compute loading at. 
	# Bring us the GRACE netfile array
	print station_list
	my_station_array     = helper_functions.read_station_file(station_list); # Structure has: name, lon, lat, starttime, endtime
	grace_data_structure = helper_functions.read_GRACE_netfile(netfile, scalefile);  # Structure has: grace_start_times,layerstrings,x_range,y_range,we,scalefactor
	return [my_station_array, grace_data_structure];



# ---------------------- COMPUTATION ------------------------ # 
def compute(output_dir, file_prefix, use_scale, disk_radius, PREM_green, threshhold_distance, my_station_array, grace_data_structure):
	# THE COMPUTATION LOOP: ONCE FOR EACH STATION
 	for i in range(len(my_station_array.station_names)):
		print "Running "+file_prefix+"computation for station "+my_station_array.station_names[i];
		output_file = output_dir+file_prefix+my_station_array.station_names[i]+"_PREM_model_ts.txt";
		compute_load_ts(my_station_array.start_time[i], my_station_array.end_time[i], my_station_array.station_lon[i], my_station_array.station_lat[i], 
			"temp_my_station.txt", "temp_GRACE_load.txt", output_file, threshhold_distance, use_scale, grace_data_structure, disk_radius, PREM_green);
	return;


def compute_load_ts(start_date, end_date, station_lon, station_lat, stationfile, loadfile, output_file, threshhold_distance, use_scale, 
	grace_data_structure, disk_radius, PREM_green):
	"""
	# For a given GPS station: 
	# For each timestamp after starting time: 
	# Write station_file
	# Write load_file
	# Compute load (mode 1 for append)
	# timestamp   --> output
	# Deformation --> output """

	for n in range(len(grace_data_structure.grace_start_times)):  # for each snapshot: 
		if grace_data_structure.grace_start_times[n]<start_date or grace_data_structure.grace_start_times[n]>end_date:
			continue;
		else:
			disk_loadfile="disk_loadfile.txt"
			helper_functions.write_load_file(grace_data_structure, n, use_scale, loadfile);  
			helper_functions.write_station_file(stationfile, station_lon, station_lat);  # assumes mode 1, just a single station. 
			helper_functions.write_time_stamp(output_file,grace_data_structure.layerstrings[n]);
			
			prem_functions.compute_prem_load(stationfile, loadfile, disk_loadfile, disk_radius, PREM_green, output_file, threshhold_distance, 1);
			# mode 1 means time series mode for a single station.
	return;



# ---------------------- OUTPUTS ------------------------ # 
def outputs(output_dir, file_prefix, my_station_array):

	# making plots
	for item in my_station_array.station_names:	
		plotroot=output_dir+file_prefix+item+"_PREM";  # where you can find/name the data
		helper_functions.make_single_plot(plotroot, item); 

	# cleaning up
	helper_functions.delete_files_matching("temp*.txt");
	helper_functions.delete_files_matching("*.pyc");
	helper_functions.delete_files_matching("disk_loadfile.txt");
	return; 




