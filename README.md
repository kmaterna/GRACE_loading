# Evaluate GRACE loads on GPS Stations
Written by Kathryn Materna, 2017

This program takes all GRACE gravity loads within a certain distance from station, 
and computes and integrates the loading effect from each cell.  
You can choose:
- to compute loads on a PREM earth structure or elastic half-space. 
- to use one of several Green's functions for loading disks on a PREM earth.
- to use which GRACE solution (JPL, GFZ, CSR, Mascons). 
- to use the scaling grid (1) or not use scaling grid (0). 

Pseudo code:  
Read configfile
Read in station information for your network: `name lon lat T1 T2`  
- For each station, 
  - For each timestep from T1 to T2,
    - load grid ---> temp input file
    - station loc -> temp input file
    - timestamp ---> station output file
    - 3D disp.  ---> station output file 
Write output textfiles with east, north, up modeled displacements
Produce output timeseries plots

### To run:  
The input string is of the format `grace_loading_driver.py config.txt`  
Make sure that haversine is on your pythonpath.  


### Example Output Plot
![CoulombCalc](https://github.com/kmaterna/GRACE_loading/blob/master/Example/RPUR_PREM_grace.png)