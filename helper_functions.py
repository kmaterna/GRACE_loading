#!/usr/env/bin python

# September 18, 2017
# Helper functions for the elastic loading time series calculations
# These are common to all calculations, no matter what Green's functions we use. 


import datetime as dt 
from scipy.io import netcdf
import collections
import os 
from subprocess import call
import glob
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 

Station_Array=collections.namedtuple('Station_Array',['station_names',
	'station_lon', 'station_lat','start_time','end_time']);
GRACE_Array=collections.namedtuple('GRACE_Array',['grace_start_times',
	'layerstrings','x_range','y_range','we','scalefactor']);



def define_io_options(network, dataset, use_scale):
	""" This simple function sets up the paths for outputs depending on what options we want (JPL vs. GFZ, scaled vs. unscaled).
		It tells the program where to find the data and where to send the results. 
	"""
	output_stem="OUTPUT/"+network;
	grace_data_stem="../GRACEDATA/DATA/"  # where the GRACE netcdf file live. 
	print "Output stem is "+output_stem+"/";

	if not os.path.exists(output_stem):
		call(['mkdir',output_stem],shell=False);
	if dataset == "JPL":
		netfile = grace_data_stem+"GRCTellus.JPL.200204_201701.LND.RL05_1.DSTvSCS1411.nc";
	if dataset == "GFZ":
		netfile = grace_data_stem+"GRCTellus.GFZ.200204_201701.LND.RL05.DSTvSCS1409.nc";
	if dataset == "CSR":
		netfile = grace_data_stem+"GRCTellus.CSR.200204_201701.LND.RL05.DSTvSCS1409.nc";
	if use_scale == 1:
		file_prefix="scaled";
	else:
		file_prefix="unscaled";	
	output_dir = output_stem+"/"+dataset+"_"+file_prefix;  # Example: "MIBB_stations/JPL_scaled/"
	if not os.path.exists(output_dir):
		call(['mkdir',output_dir],shell=False);
	scalefile="../GRACEDATA/DATA/CLM4.SCALE_FACTOR.DS.G300KM.RL05.DSTvSCS1409.nc";

	return [netfile, scalefile, output_dir+"/", file_prefix+"_"];



def read_GRACE_netfile(netfile, scalefile):
	"""
	Meant for reading a netcdf file that has many time-slices. 
	"""

	f = netcdf.netcdf_file(netfile,'r');
	layernames = f.input_filename;
	layernames = layernames.split('\r');
	layerstrings = titlestrings(layernames);  # Format: 01-Apr-2002_30-Jun-2002
	grace_start_times = beginepochs(layernames);    # Format: 2014.254000

	# Get the x_range and y_range (both just two numbers, like [39.38 39.86] degrees N)
	x_range = f.variables['lon'];
	x_range = x_range[:].copy();
	y_range = f.variables['lat'];
	y_range = y_range[:].copy();
	z_range = f.variables['time'];
	z_range = z_range[:].copy();

	# Get the grids and labels for gravity / water thickness. 
	we = f.variables['lwe_thickness'];
	we = we[:].copy();

	# Get the GRACE scale factor corrections. 
	sf = netcdf.netcdf_file(scalefile,'r');
	scalefactor=sf.variables['SCALE_FACTOR'];   # size = 360x180; time-independent. 
	scalefactor=scalefactor[:].copy();

	grace_data_structure=GRACE_Array(grace_start_times=grace_start_times, layerstrings=layerstrings, 
		x_range=x_range, y_range=y_range, we=we, scalefactor=scalefactor);

	return grace_data_structure; 



def read_station_file(station_list):

	ifile=open(station_list,'r');
	station_names=[]; station_lon=[]; station_lat=[]; start_time=[]; end_time=[];
	grace_begin=2002.000;
	grace_end=2017.1;
	grace_default_start=2012.000;
	for line in ifile:
		temp=line.split();
		station_names.append(temp[0]);
		station_lon.append(temp[1]);
		station_lat.append(temp[2]);
	
		# Decide on the time range for the computation at each GPS station. 
		# If time constraints are given: 
		if len(temp)>4:
			gps_start_time=float(temp[3]);
			gps_end_time  =float(temp[4]);
			if gps_start_time<grace_begin:
				gps_start_time=grace_begin;
				# We can't compute GRACE before 2002.000, so ignore all data before that. 
			gps_start_time=max(gps_start_time-1.5, grace_begin);
			gps_end_time  =min(gps_end_time+1.5,grace_end);
			# Here we adopt the slightly padded range of the GPS data for our GRACE calculation

		else: 		# If time constraints are not given, replace them with default time window. 
			gps_start_time=grace_default_start;
			gps_end_time=grace_end;

		start_time.append(gps_start_time);
		end_time.append(gps_end_time);
		my_station_array=Station_Array(station_names=station_names,station_lon=station_lon, 
			station_lat=station_lat,start_time=start_time,end_time=end_time);
	return my_station_array; 


def write_load_file(grace_data_structure, n, use_scale,loadfile):
	"""
	This reads the n'th slice from the netfile [X x Y x n], multiplies by the scale file, and writes the load into a load file. 
	If use_scale is set, then we multiply the GRACE timeslice by the scale factor grid. 
	This is a global load file. 
	"""	
	# Multiply by scale factor (optional).  Write out to file.  Format = lon, lat, xspacing, yspacing, load
	outfile=open(loadfile,'w');
	
	snapshot=grace_data_structure.we[n,:,:];   # Get a snapshot
	s1=np.shape(snapshot);   # fill in NaN's as necessary over the water mask. 
	mult_scalefactor=grace_data_structure.scalefactor;
	if use_scale==0:  # If no scaling grid applied: 
		mult_scalefactor=np.ones(s1);  # if we are not using scale factor, we just multiply everything by 1. 

	for i in range(s1[0]):
		for j in range(s1[1]):
			if snapshot[i,j]<5000:        # Write the numeric values for pixels not over water. 
				pressure = snapshot[i,j]*mult_scalefactor[i,j]*(1.0/100)*9.81*1000.0;  
				# rho*g*h (h converted cm --> meters)   
				outfile.write("%f %f 0.5 0.5 %f\n" %(grace_data_structure.x_range[j], grace_data_structure.y_range[i], pressure));  
				# write an ascii file for GMT.
				# The 0.5 and 0.5 are for the GRACE 1-degree-by-1-degree grids. 
	outfile.close();
	return;



def write_station_file(stationfile, station_lon, station_lat):
	ifile=open(stationfile,'w');
	ifile.write("%s %s 0\n" %(station_lon, station_lat));
	ifile.close();
	return;

def write_time_stamp(outfile,timestamp):
	ifile=open(outfile,'a');
	ifile.write(str(timestamp)+' ');
	ifile.close();
	print timestamp
	return;

def get_decyear(yyyyddd):
	myday=dt.datetime.strptime(yyyyddd,"%Y%j");
	year=float(yyyyddd[0:4]);
	jday=float(myday.strftime("%j"));
	decyear=year+jday/365.24;
	return decyear;

def get_ddmmyyyy(yyyyddd):
	myday=dt.datetime.strptime(yyyyddd,"%Y%j");
	year=yyyyddd[0:4];
	month=myday.strftime("%b");
	dd=myday.strftime("%d");
	formatstring=dd+"-"+month+"-"+year;
	return formatstring;

def get_datetimes(tsfile):
	ifile=open(tsfile,'r');
	dateobjects=[];
	for line in ifile:
		temp=line.split();
		raw_string=temp[0];   # This is in the format 01-Jan-2012_31-Jan-2012
		datestring=raw_string[0:11];
		myobject=dt.datetime.strptime(datestring,'%d-%b-%Y');
		dateobjects.append(myobject);
	return dateobjects;

def titlestrings(layernames):
	"""
	Take a set of names like './GSM/GSM-2_2014121-2014151_0031_JPLEM_0001_0005\n'
	and convert them into '01-Apr-2014_30-Jun-2014'
	"""
	layerstrings=[];
	for item in layernames:
		yyyyddd1=item[12:19];  # start date
		yyyyddd2=item[20:27];  # end date
		converted_name1=get_ddmmyyyy(yyyyddd1);
		converted_name2=get_ddmmyyyy(yyyyddd2);
		layerstrings.append(converted_name1+'_'+converted_name2);
	return layerstrings;

def beginepochs(layernames):
	"""
	Take a set of names like './GSM/GSM-2_2014121-2014151_0031_JPLEM_0001_0005\n'
	and convert them into '2014.2500'
	"""
	start_times=[];
	for item in layernames:
		yyyyddd1=item[12:19];  # start date
		start_times.append(get_decyear(yyyyddd1));
	return start_times;

def delete_files_matching(match_string):
	clear_list=glob.glob(match_string);
	for item in clear_list:
		call(['rm',item],shell=False); 	
	return;

def make_single_plot(plotroot, station_name):
	tsfile=plotroot+'_model_ts.txt';
	[x, y, z, u, v, w] = np.loadtxt(tsfile,usecols=range(1,7),unpack=True);
	datetimes=get_datetimes(tsfile);
	t = matplotlib.dates.date2num(datetimes);
	plt.figure();
	plt.plot_date(t,w,'b--');
	plt.plot_date(t,u,'r');
	plt.plot_date(t,v,'g');
	plt.plot_date(t,w,'b.',markersize=15);
	plt.ylabel('Displacement (mm)');
	plt.title(station_name+' Model Time Series');
	plt.legend(['Vertical','East','North']);
	plt.savefig(plotroot+'.jpg');
	plt.close();
	return;


