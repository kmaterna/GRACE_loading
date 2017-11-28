#!/usr/env/bin python

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

