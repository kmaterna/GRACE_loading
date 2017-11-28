Evaluate GRACE loads on GPS Stations
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
