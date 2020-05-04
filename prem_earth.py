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
import os, sys
import matplotlib.pyplot as plt 
import prem_functions
import helper_functions

# Major interface: 
# output = mass_loading_python_implementation(GRACE_array, stations_array);

def prem_earth_grace_timeseries(params):
	netfile, scalefile = configure_filepaths(params);
	station_obj     = helper_functions.read_station_file(params.input_file); # named tuple
	grace_data_structure = helper_functions.read_GRACE_netfile(netfile, scalefile);  # named tuple
	
	computer = pyimpl_factory(params);
	for i in range(len(station_obj.name)):  # once for each station
		single_station=helper_functions.Station_Array(name=[station_obj.name[i]],lon=[station_obj.lon[i]],lat=[station_obj.lat[i]],start_time=[station_obj.start_time[i]],end_time=[station_obj.end_time[i]]);
		dates, disp_x, disp_y, disp_v = computer(grace_data_structure, single_station);
		outputs(dates, disp_x, disp_y, disp_v, params, single_station);
	return;


def pyimpl_factory(params):
	# The implementation of our calculation in python. 
	# This can be different for the half-space and C case. 
	GF = helper_functions.read_spherical_greensfunc(params.greensfunc);  
	disk_radius=6.206;

	def mass_loading_TS_python_implementation(GRACE_array, single_station):
		dates_used=[]; disp_x=[]; disp_y=[]; disp_v=[];
		for n in range(len(GRACE_array.grace_central_times)):  # for each snapshot: 
			if GRACE_array.grace_central_times[n]<single_station.start_time[0] or GRACE_array.grace_central_times[n]>single_station.end_time[0]:
				continue;
			else:
				print("Calculating %s for station %s" % (GRACE_array.grace_central_times[n],single_station.name[0]) );
				GRACE_single_loads = helper_functions.package_loads_from_GRACE(GRACE_array, n, params.scaling);  # more basic structure of the loads
				x, y, z = prem_functions.compute_prem_load(single_station, GRACE_single_loads, disk_radius, GF, params.max_distance); 
				disp_x.append(x[0]); 
				disp_y.append(y[0]); 
				disp_v.append(z[0]);  # single station must be unpacked. 
				dates_used.append(GRACE_array.grace_central_times[n]);
		return dates_used, disp_x, disp_y, disp_v;

	return mass_loading_TS_python_implementation;

# ---------------------- CONFIGURE ------------------------ # 
def configure_filepaths(params):
	# Sending the right datafile into the computation. 
	if params.gracetype=="Mascon":
		netfile = params.mascon_datafile;
		scalefile = params.mascon_scalefile;
	elif params.gracetype=="TELLUS":
		if params.datasource=="JPL":
			netfile = params.tellus_jpl_datafile;
		elif params.datasource=="GFZ":
			netfile = params.tellus_gfz_datafile;
		elif params.datasource=="CSR":
			netfile = params.tellus_csr_datafile;
		else:
			print("ERROR! Your TELLUS datasource is not valid [JPL, CSR, GFZ]. Please try again");
			sys.exit(0);
		scalefile=params.tellus_scalefile;
	return netfile, scalefile;


# ---------------------- OUTPUTS ------------------------ # 
def outputs(dates, disp_x, disp_y, disp_v, params, single_station):
	call(['mkdir','-p',params.output_dir],shell=False);
	output_dir = params.output_dir+"/"+params.datasource+"_"+params.gracetype+"_"+str(params.scaling)+"/";  # Example: "MIBB/JPL_TELLUS_0/"
	call(['mkdir','-p',output_dir],shell=False);
	call(['cp','config.txt',output_dir],shell=False);  # in case you used config.txt as your config file. 
	output_file = output_dir+single_station.name[0]+"_PREM_model_ts.txt";
	ofile=open(output_file,'w');
	for i in range(len(dates)):
		ofile.write(dt.datetime.strftime(dates[i],"%Y-%m-%d"));
		ofile.write(" %f %f %.5f %.5f %.5f\n" % (single_station.lon[0], single_station.lat[0], disp_x[i], disp_y[i], disp_v[i]) );
	ofile.close();
	# # making plots
	plt.figure()
	plt.plot_date(dates,disp_v,'b--');
	plt.plot_date(dates,disp_x,'r');
	plt.plot_date(dates,disp_y,'g');
	plt.plot_date(dates,disp_v,'b.',markersize=15);
	plt.ylabel('Displacement (mm)');
	plt.legend(['Vertical','East','North']);	
	plt.savefig(output_dir+single_station.name[0]+"_PREM_grace.png");
	return; 

def temp_plot(mydate):
	[lon, lat, load] = np.loadtxt("temp_global_GRACE_load.txt",unpack=True,usecols=(0,1,4));
	plt.figure(dpi=300);
	plt.scatter(lon, lat, c=load, s=0.5, vmin=-3000, vmax=3000, cmap="jet");
	plt.colorbar();
	plt.savefig(dt.datetime.strftime(mydate,"%Y%m%d")+"_loads.png");


