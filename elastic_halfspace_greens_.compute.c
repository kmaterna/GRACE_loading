#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "elastic_halfspace_load.h"

// Compute the deformation of an elastic halfspace from a set of surface load rectangles!
// This code expects:
// 1. an input file with the description of the load rectangles (lon, lat, a(half-width, degrees), b(half-degrees), pressure). 
// 2. an input file with the list of observations stations (lon, lat, depth) . 
// 3. Will make an output file with the u, v, w at each observation station. 
//
// Psuedocode: 
// Parse inputs. 
// For each line in obs_station_file:
//   For each load cell:
//     Convert position into centered cartesian meters
//     Compute u_disp, v_disp, w_disp if the load is close enough (less than threshhold_distace)
//     Write outputs. 
// Compile with this line: 
// gcc -o exec_name source_name.c 
// Call this code in the following manner:
// ./exec_name stationfile loadfile outputfile lame mu threshhold_distance ts_or_image
// ts_or_image = 1: time series (output file append-mode)
// ts_or_image = 2: image (write output file once)

// Function prototypes
double get_x_meters(double load_lon, double load_lat, double sta_lon, double sta_lat);
double get_y_meters(double load_lon, double load_lat, double sta_lon, double sta_lat);
double deg2rad(double deg);
double rad2deg(double rad);
double distance_km(double origin_lon, double origin_lat, double dest_lon, double dest_lat);
double calculate_initial_compass_bearing(double origin_lon, double origin_lat, double dest_lon, double dest_lat);

int main(int argc, char *argv[]){

	char stationfile[90];
	char loadfile[90];
	char outputfile[90];
	double lame, mu, threshhold_distance;
	double ts_or_image;
	double u_displacement, v_displacement, w_displacement; 
	double sta_lon, sta_lat, z; 
	double x, y;
	double load_lon, load_lat, adeg, bdeg, p, ameters, bmeters;
	char csta_lon[20], csta_lat[20], csta_z[20]; 
	char cload_lon[20], cload_lat[20], cadeg[20], cbdeg[20], char_p[20]; 
	char * temp;

	char * stat_line = NULL;
	char * load_line = NULL;
    size_t len = 0;
    size_t lenload = 0;
    ssize_t read, readload;
	FILE * statfile;
	FILE * lfile; 
	FILE * ofile; 
	
	// Parse and cleanup the inputs. 
	if (argc < 8) {
		printf("ERROR: Not the right number of inputs to this function. \n");
		printf("Please call this with : stationfile loadfile outputfile lame mu threshhold_distance ts_or_image \n");
		exit(1);		
	}
	strcpy(stationfile,argv[1]);
	strcpy(loadfile,argv[2]);
	strcpy(outputfile,argv[3]);
	lame = strtod(argv[4],NULL);
	mu = strtod(argv[5],NULL);
	threshhold_distance = strtod(argv[6],NULL);
	ts_or_image = strtod(argv[7],NULL);
	int mode_flag = (int) ts_or_image;

    // Open the output file in time series mode
    if (mode_flag == 1) {
    	ofile = fopen(outputfile,"a");
    	printf("Opening new output file in append mode for time series...\n");
	}
	
	// Open the output file in image mode. 
    else {
    	ofile = fopen(outputfile,"w");
    	printf("Opening output file in write mode for single image...\n");
	}

	// // Begin to read the input files. 
	statfile = fopen(stationfile,"r");
    if (statfile == NULL) {
        exit(EXIT_FAILURE);
        printf("ERROR: Station file was not opened.\n");
    }

	while ((read = getline(&stat_line, &len, statfile)) != -1){  // for each station we want to see... 
		// The program always gets to here. 
		temp=strtok(stat_line, " ");
		strcpy(csta_lon,temp);
		sta_lon = strtod(csta_lon, NULL);
		temp=strtok(NULL," ");
		strcpy(csta_lat, temp);
		sta_lat = strtod(csta_lat, NULL);
		temp=strtok(NULL," ");
		strcpy(csta_z,temp);
		z   = strtod(csta_z,   NULL);
		// printf("Lon Lat Z are: %f %f %f \n", sta_lon, sta_lat, z);
		// Now we have a station latitude and longitude. 

		// Initialize variables
		u_displacement = 0;
		v_displacement = 0;
		w_displacement = 0;

		// Integrate over load cells. 
		lfile = fopen(loadfile,"r");
		if (lfile == NULL) {
		    exit(EXIT_FAILURE);
		    printf("ERROR: Load file was not opened.\n");
		}
		// for each cell in the applied load
		while ((readload = getline(&load_line, &lenload, lfile)) != -1){  // for each cell in the applied load
			temp=strtok(load_line, " ");
			strcpy(cload_lon,temp);
			load_lon = strtod(cload_lon, NULL);
			temp=strtok(NULL," ");
			strcpy(cload_lat, temp);
			load_lat = strtod(cload_lat, NULL);
			temp=strtok(NULL," ");
			strcpy(cadeg,temp);
			adeg  = strtod(cadeg,   NULL);
			temp=strtok(NULL," ");
			strcpy(cbdeg,temp);
			bdeg   = strtod(cbdeg,   NULL);
			temp=strtok(NULL," ");
			strcpy(char_p,temp);  
			p = strtod(char_p,   NULL);
			// printf("Lon Lat a b p are: %f %f %f %f %f\n", load_lon, load_lat, adeg, bdeg, p);

			// a and b in meters
			ameters = (adeg/1.0)*111000*cos(pi*sta_lat/180);     // half-width degrees to meters, east-west rectangle. 
			bmeters = (bdeg/1.0)*111000;  // half-width degrees to meters, north-south rectangle. 

			x = get_x_meters(load_lon, load_lat, sta_lon, sta_lat);
			// printf("%f %f %f %f \n",load_lon, load_lat, sta_lon, sta_lat);
			// printf("distance: %f\n",x);
			y = get_y_meters(load_lon, load_lat, sta_lon, sta_lat);

			if ( sqrt( pow(x,2) + pow(y,2)) < threshhold_distance ){
				u_displacement+=1000.0 * u(x, y, z, lame, mu, p, ameters, bmeters);
				v_displacement+=1000.0 * v(x, y, z, lame, mu, p, ameters, bmeters);
				w_displacement+=1000.0 * w(x, y, z, lame, mu, p, ameters, bmeters);
			}

		}
		fclose(lfile);  // done with integrating cells of the applied load. 
		// write to the output file. 
		fprintf(ofile,"%f %f %f %f %f %f\n", sta_lon, sta_lat, z, u_displacement, v_displacement, w_displacement);
		// Now go back and compute everything for a new station. 

	} // done calculating and writing out for all desired stations.
	fclose(statfile);
	fclose(ofile);
	return(0);
}



double distance_km(double origin_lon, double origin_lat, double dest_lon, double dest_lat){
    /*
    Computes the distance [in km] between origin [lat1, lon1] and destination [lat2, lon2]. 
    Lat-lon must be specified in that order. 
    Haversine Formula.
    */
    double radius, dlon, dlat, a, c, d;
    radius = 6371; // km, average radius of earth. 
    dlat = deg2rad(dest_lat-origin_lat);
    dlon = deg2rad(dest_lon-origin_lon);
    a = sin(dlat/2) * sin(dlat/2) + cos(deg2rad(origin_lat)) * cos(deg2rad(dest_lat)) * sin(dlon/2) * sin(dlon/2);
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    d = radius * c;
    return d;
}


double calculate_initial_compass_bearing(double origin_lon, double origin_lat, double dest_lon, double dest_lat){
    /* """
    Calculates the bearing between two points.
    The formulae used is the following:
        theta = atan2(sin(delta_long).cos(lat2),cos(lat1).sin(lat2) - sin(lat1).cos(lat2).cos(delta_long))
    :Parameters:
      - `pointA: The tuple representing the latitude/longitude for the
        first point. Latitude and longitude must be in decimal degrees
      - `pointB: The tuple representing the latitude/longitude for the
        second point. Latitude and longitude must be in decimal degrees
    :Returns:
      The bearing in degrees (CW from north, just like strike)
    :Returns Type: float
    """ */

    double lat1, lat2, diffLong, x, y, initial_bearing, compass_bearing;
    lat1 = deg2rad(origin_lat);
    lat2 = deg2rad(dest_lat);

    diffLong = deg2rad(dest_lon - origin_lon);

    x = sin(diffLong) * cos(lat2); 
    y = cos(lat1) * sin(lat2) - (sin(lat1) * cos(lat2) * cos(diffLong));
    initial_bearing = atan2(x, y);

    // Now we have the initial bearing but math.atan2 return values
    // from -180 to + 180 which is not what we want for a compass bearing
    // The solution is to normalize the initial bearing as shown below
    initial_bearing = rad2deg(initial_bearing);
    compass_bearing = fmod((initial_bearing + 360.0), 360.0);
    return compass_bearing;
}

double get_x_meters(double load_lon, double load_lat, double sta_lon, double sta_lat){
    /* """
    Distance between two latitude/longitude pairs, 
    given in x-distance and y-distance in meters
    (assuming flat surface between the points)
    """ */
    double radius, bearing, azimuth, x;
    radius = distance_km(load_lon, load_lat, sta_lon, sta_lat);
    bearing = calculate_initial_compass_bearing(load_lon, load_lat, sta_lon, sta_lat);
    azimuth = 90 - bearing;
    x = radius * cos(deg2rad(azimuth)) * 1000;
    return x;
}

double get_y_meters(double load_lon, double load_lat, double sta_lon, double sta_lat){
    /*"""
    Distance between two latitude/longitude pairs, 
    given in x-distance and y-distance in meters
    (assuming flat surface between the points)
    """ */
    double radius, bearing, azimuth, y;
    radius = distance_km(load_lon, load_lat, sta_lon, sta_lat);
    bearing = calculate_initial_compass_bearing(load_lon, load_lat, sta_lon, sta_lat);
    azimuth = 90 - bearing;
    y = radius * sin(deg2rad(azimuth)) * 1000;
    return y;
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts decimal degrees to radians             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double deg2rad(double deg) {
  return (deg * pi / 180);
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts radians to decimal degrees             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double rad2deg(double rad) {
  return (rad * 180 / pi);
}
