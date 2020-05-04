#!/usr/env/bin python

import numpy as np 
from subprocess import call 
import datetime as dt 
import os
import sys
import matplotlib
import matplotlib.pyplot as plt 
import helper_functions
import glob




def plot_many_GRACE_types(network):
	[station_names, output_dir, model] = configure(network);
	outputs(output_dir, station_names);
	return;

# ---------------------- CONFIGURE ------------------------ # 
def configure(network):
	""" This simple function sets up the paths for outputs depending on what options we want (JPL vs. GFZ, scaled vs. unscaled).
		It tells the program where to find the data and where to send the results. 
	"""
	input_file = "INPUT/"+network+".txt"
	station_names=[];
	ifile=open(input_file);
	for line in ifile:
		station_names.append(line.split()[0])

	output_dir="OUTPUT/"+network+"/";
	print "Output stem is "+output_dir;
	if not os.path.exists(output_dir):
		print "Problem! No Output Directory!"

	model="PREM"
		
	return [ station_names, output_dir , model];


# ---------------------- OUTPUTS ------------------------ # 
def outputs(output_dir, station_names, model):
	for item in station_names:	
		filelist=glob.glob(output_dir+"*/*"+item+"_"+model+"_model_ts.txt");
		print "Existing files are...";
		print filelist;
		make_multi_GRACE_plot(filelist, item,output_dir);
	return;

def get_label_from_tsfile(tsfile):
	if "GFZ" in tsfile:
		procenter="GFZ"
	elif "JPL" in tsfile:
		procenter="JPL"
	elif "CSR" in tsfile:
		procenter="CSR"
	if "_scaled" in tsfile:
		scaling="scaled";
	if "_unscaled" in tsfile:
		scaling="unscaled"
	return procenter+" "+scaling;


def make_multi_GRACE_plot(filelist, station_name, output_dir):

	color_array=['blue','skyblue','red','indianred','green','lightgreen'];
	plt.figure();
	
	for i in range(len(filelist)):
		tsfile=filelist[i];
		mylabel=get_label_from_tsfile(tsfile);
		[x, y, z, u, v, w] = np.loadtxt(tsfile,usecols=range(1,7),unpack=True);
		datetimes=helper_functions.get_datetimes(tsfile);
		t = matplotlib.dates.date2num(datetimes);
		plt.plot_date(t,w,marker='.',linestyle='--',color=color_array[i],markersize=15,label=mylabel);

	plt.ylabel('Displacement (mm)');
	plt.title(station_name+' Model Time Series');
 	plt.legend();
	plt.savefig(output_dir+station_name+'.jpg');
	plt.close();
	return;


if __name__=="__main__":
	# OPTIONS
	# The options that you might want to change. 	
	if len(sys.argv)==2:
		network=str(sys.argv[1]);
	else:
		print "You should provide manual arguments to this script\n such as 'MIBB'\n"; exit(1);

	plot_many_GRACE_types(network);

