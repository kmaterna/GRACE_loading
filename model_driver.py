#!/usr/env/bin python
"""
EVALUATE GRACE LOADS ON GPS STATIONS
Written by Kathryn Materna, 2017
This program takes all GRACE gravity loads within a certain distance from station, 
and computes the loading effect from each 1-degree-by-1-degree cell. 
You can choose:
	-to compute loads on a PREM earth structure or elastic half-space. 
	-to use which GRACE solution (JPL, GFZ, CSR, Mascons). 
	-to use the scaling grid (1) or not use scaling grid (0). 
Pseudo code: 
Read in station information for your network: name, lon, lat, T1, T2. 
For each station, 
 For each timestep from T1 to T2,
   compute 3D loading displacement
To run: "python model_driver.py config.txt"
Make sure that the green's functions and haversine are on your pythonpath. 
"""

import sys
import argparse, configparser
import collections
import helper_functions
import prem_earth

Params = collections.namedtuple("Params",["network","input_file","output_dir","grace_dir",
	"mascon_datafile","mascon_scalefile","tellus_jpl_datafile","tellus_csr_datafile","tellus_gfz_datafile",
	"tellus_scalefile","grace_dates_file","gracetype","max_distance","datasource","greensfunc","scaling"]);

def welcome_and_parse():
	print("\n\nWelcome to a simple forward modeling tool for calculating GRACE loading at GPS points. ");
	parser = argparse.ArgumentParser(description='Run GRACE load models in Python', epilog='\U0001f600 \U0001f600 \U0001f600 ');
	parser.add_argument('config',type=str,help='name of config file for calculation. Required.')
	args = parser.parse_args()
	print("Config file:",args.config);
	return args;

def configure_calc(config_file):
	configobj=configparser.ConfigParser();
	configobj.optionxform = str # make the config file case-sensitive
	configobj.read(config_file);

	# IO parameters
	network = configobj.get('io-config','network')
	input_file=configobj.get('io-config','gps_infile');
	output_dir=configobj.get('io-config','output_dir');
	grace_dir =configobj.get('io-config','GRACE_data_dir');

	mascon_datafile = grace_dir+configobj.get('io-config','mascon_datafile');
	mascon_scalefile = grace_dir+configobj.get('io-config','mascon_scalefile');
	tellus_jpl_datafile = grace_dir+configobj.get('io-config','tellus_jpl_datafile');
	tellus_csr_datafile = grace_dir+configobj.get('io-config','tellus_csr_datafile');
	tellus_gfz_datafile = grace_dir+configobj.get('io-config','tellus_gfz_datafile');
	tellus_scalefile = grace_dir+configobj.get('io-config','tellus_scalefile');
	grace_dates_file = grace_dir+configobj.get('io-config','grace_dates_file');

	# Computation parameters
	greensfunc = configobj.get('calc-config','greensfunc');
	gracetype = configobj.get('calc-config','gracetype');
	max_distance = configobj.getfloat('calc-config','max_distance');
	datasource = configobj.get('calc-config','datasource');	
	scaling = configobj.getint('calc-config','scaling');

	MyParams = Params(network=network, input_file=input_file, output_dir=output_dir, grace_dir=grace_dir,
		mascon_datafile=mascon_datafile, mascon_scalefile=mascon_scalefile, tellus_jpl_datafile=tellus_jpl_datafile, 
		tellus_csr_datafile=tellus_csr_datafile, tellus_gfz_datafile=tellus_gfz_datafile, tellus_scalefile=tellus_scalefile, 
		grace_dates_file=grace_dates_file, gracetype=gracetype,
		max_distance=max_distance, datasource=datasource, greensfunc=greensfunc, scaling=scaling);
	print(MyParams);
	print("\n");
	return MyParams;	


if __name__=="__main__":

	# The options that you might want to change. 
	args = welcome_and_parse();
	params = configure_calc(args.config);
	prem_earth.prem_earth_grace_timeseries(params);


