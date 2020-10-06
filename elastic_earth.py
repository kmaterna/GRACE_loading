#!/usr/env/bin python

from subprocess import call
import helper_functions


def elastic_earth_computation(network, dataset, use_scale):
    [netfile, scalefile, station_file, output_dir, file_prefix, mu, lame, threshhold_distance] = configure(network,
                                                                                                           dataset,
                                                                                                           use_scale);
    [my_station_array, grace_data_structure] = get_inputs(station_file, netfile, scalefile);
    compute(output_dir, file_prefix, use_scale, mu, lame, threshhold_distance, my_station_array, grace_data_structure);
    outputs(output_dir, file_prefix, my_station_array);
    return;


# ---------------------- CONFIGURE ------------------------ # 
def configure(network, dataset, use_scale):
    # GET STARTED WITH SET-UP AND CLEAN-UP
    input_file = "INPUT/" + network + ".txt"
    [netfile, scalefile, output_dir, file_prefix] = helper_functions.define_io_options(network, dataset, use_scale);

    helper_functions.delete_files_matching(
        output_dir + file_prefix + '*_halfspace_model_ts.txt');  # get started with a clean directory.

    call('gcc -o a elastic_halfspace_greens_compute.c', shell=True);  # compile C code

    # Other parameters I don't often change.
    mu = 30 * 1000 * 1000 * 1000;  # assuming poisson's ratio is 0.25, then mu = lame elastic parameter
    lame = mu;
    threshhold_distance = 2000 * 1000;  # in m. Further than 2000km away, we assume the elastic loads don't matter.
    return [netfile, scalefile, input_file, output_dir, file_prefix, mu, lame, threshhold_distance];


# ---------------------- INPUTS ------------------------ # 
def get_inputs(station_list, netfile, scalefile):
    # Bring us arrays of lat, lon, and time ranges that we need to compute loading at.
    # Bring us the GRACE netfile array
    print(station_list)
    my_station_array = helper_functions.read_station_file(
        station_list);  # Structure has: name, lon, lat, starttime, endtime
    grace_data_structure = helper_functions.read_GRACE_netfile(netfile, scalefile)
    return [my_station_array, grace_data_structure];


# ---------------------- COMPUTATION ------------------------ #
def compute(output_dir, file_prefix, use_scale, mu, lame, threshhold_distance, my_station_array, grace_data_structure):
    # THE COMPUTATION LOOP
    for i in range(len(my_station_array.station_names)):
        print("Running " + file_prefix + "computation for station " + my_station_array.station_names[i]);
        output_file = output_dir + file_prefix + my_station_array.station_names[i] + "_halfspace_model_ts.txt";
        compute_load_ts(my_station_array.start_time[i], my_station_array.end_time[i], my_station_array.station_lon[i],
                        my_station_array.station_lat[i],
                        "temp_my_station.txt", "temp_GRACE_load.txt", output_file, mu, lame, threshhold_distance,
                        use_scale, grace_data_structure);
    return;


def compute_load_ts(start_date, end_date, station_lon, station_lat, stationfile, loadfile, output_file, lame, mu,
                    threshhold_distance, use_scale, grace_data_structure):
    """
    # For a given GPS station:
    # For each timestamp after starting time:
    # Write station_file
    # Write load file
    # Compute load (mode 1 for append)
    # timestamp   --> output
    # Deformation --> output """

    for n in range(len(grace_data_structure.grace_start_times)):  # for each snapshot:
        if grace_data_structure.grace_start_times[n] < start_date or grace_data_structure.grace_start_times[n] > end_date:
            continue;
        else:
            helper_functions.write_load_file(grace_data_structure, n, use_scale, loadfile);
            helper_functions.write_station_file(stationfile, station_lon,
                                                station_lat);  # write a file with one station's coordinate.
            helper_functions.write_time_stamp(output_file, grace_data_structure.layerstrings[
                n]);  # write the time stamp into the output file.

            # Do the computation: Mode 1 = time series (append to the output file).
            # Mode 2 = image (just write the output file).
            call(['./a', stationfile, loadfile, output_file, str(lame), str(mu), str(threshhold_distance), str(1)],
                 shell=False);
    return;


# ---------------------- OUTPUTS ------------------------ #
def outputs(output_dir, file_prefix, my_station_array):
    for item in my_station_array.station_names:
        plotroot = output_dir + file_prefix + item + "_halfspace";  # where you can find/name the data
        helper_functions.make_single_plot(plotroot, item);

    helper_functions.delete_files_matching("temp*.txt");
    helper_functions.delete_files_matching("*.pyc");
    helper_functions.delete_files_matching("a");
    return;
