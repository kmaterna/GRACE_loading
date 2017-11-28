#!/usr/env/bin python
"""
EVALUATE GRACE LOADS ON GPS STATIONS
Written by Kathryn Materna, 2017
This program takes all GRACE gravity loads within a certain distance from station, 
and computes the loading effect from each 1-degree-by-1-degree cell. 
You can choose:
	-to compute loads on a PREM earth structure or elastic half-space. 
	-to use which GRACE solution (JPL, GFZ, CSR). 
	-to use the scaling grid (1) or not use scaling grid (0). 

Pseudo code: 
Go into INPUT/my_network.txt. 
Read in station information for your network: name, lon, lat, T1, T2. 
For each station, 
 For each timestep from T1 to T2,
   load grid ---> temp input file
   station loc -> temp input file
   timestamp ---> station output file
   3D disp.  ---> station output file 
To run: 
The input string is "python model_driver.py MIBB JPL 0"
Make sure that the green's functions and haversine are on your pythonpath. 
"""

import sys
import elastic_earth
import prem_earth

if __name__=="__main__":

	# The options that you might want to change. 	
	if len(sys.argv)==4:
		network=str(sys.argv[1]);
		dataset=str(sys.argv[2]);  # options are JPL or GFZ of CSR. 
		use_scale=int(sys.argv[3]);  # scale factor is something that we use to un-do the automatic Gaussian filtering in GRACE data. 	
	else:  # if you're controlling the script manually:
		print "Error! please configure run-string like: 'python model_driver.py MIBB JPL 0' \n" 
		exit(0)


	# elastic_earth.elastic_earth_computation(network, dataset, use_scale);
	prem_earth.prem_earth_grace_timeseries(network, dataset, use_scale);

