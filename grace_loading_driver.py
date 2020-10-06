#!/usr/bin/env python
"""
EVALUATE GRACE LOADS ON GPS STATIONS
Written by Kathryn Materna, 2017
This program takes all GRACE gravity loads (MASCON or TELLUS) within a certain distance from station, 
and computes the loading effect from each cell. 
You can choose:
    -to compute loads on a PREM earth structure or elastic half-space.
    -to use which GRACE solution (JPL, GFZ, CSR, Mascons).
    -to use the scaling grid (1) or not use scaling grid (0).
Inputs: 
Read in station information for your network: name, lon, lat, T1, T2. 
For each station,
    For each timestep from T1 to T2,
    compute 3D loading displacement
To run: "grace_loading_driver.py config.txt"
Make sure that haversine is on your pythonpath. 
"""

import argparse
import parse_configfile
import prem_earth


def welcome_and_parse():
    print("\n\nWelcome to a forward modeling tool for calculating GRACE loading at GPS points. ");
    parser = argparse.ArgumentParser(description='Run GRACE load models in Python',
                                     epilog='\U0001f600 \U0001f600 \U0001f600 ');
    parser.add_argument('configfile', type=str, help='name of config file for calculation. Required.')
    args = parser.parse_args()
    print("Config file:", args.configfile);
    return args;


if __name__ == "__main__":
    # The options that you might want to change.
    args = welcome_and_parse();
    params = parse_configfile.configure_calc(args.configfile);
    prem_earth.prem_earth_grace_timeseries(params);
