# Read the config file

import configparser
import collections

Params = collections.namedtuple("Params", ["network", "input_file", "output_dir", "grace_dir",
                                           "mascon_datafile", "mascon_scalefile", "tellus_jpl_datafile",
                                           "tellus_csr_datafile", "tellus_gfz_datafile",
                                           "tellus_scalefile", "grace_dates_file", "gracetype", "max_distance",
                                           "datasource", "greensfunc", "scaling"]);


def configure_calc(config_file):
    configobj = configparser.ConfigParser();
    configobj.optionxform = str  # make the config file case-sensitive
    configobj.read(config_file);

    # IO parameters
    network = configobj.get('io-config', 'network')
    input_file = configobj.get('io-config', 'gps_infile');
    output_dir = configobj.get('io-config', 'output_dir');
    grace_dir = configobj.get('io-config', 'GRACE_data_dir');

    mascon_datafile = grace_dir + configobj.get('io-config', 'mascon_datafile');
    mascon_scalefile = grace_dir + configobj.get('io-config', 'mascon_scalefile');
    tellus_jpl_datafile = grace_dir + configobj.get('io-config', 'tellus_jpl_datafile');
    tellus_csr_datafile = grace_dir + configobj.get('io-config', 'tellus_csr_datafile');
    tellus_gfz_datafile = grace_dir + configobj.get('io-config', 'tellus_gfz_datafile');
    tellus_scalefile = grace_dir + configobj.get('io-config', 'tellus_scalefile');
    grace_dates_file = grace_dir + configobj.get('io-config', 'grace_dates_file');

    # Computation parameters
    greensfunc = configobj.get('calc-config', 'greensfunc');
    gracetype = configobj.get('calc-config', 'gracetype');
    max_distance = configobj.getfloat('calc-config', 'max_distance');
    datasource = configobj.get('calc-config', 'datasource');
    scaling = configobj.getint('calc-config', 'scaling');

    MyParams = Params(network=network, input_file=input_file, output_dir=output_dir, grace_dir=grace_dir,
                      mascon_datafile=mascon_datafile, mascon_scalefile=mascon_scalefile,
                      tellus_jpl_datafile=tellus_jpl_datafile,
                      tellus_csr_datafile=tellus_csr_datafile, tellus_gfz_datafile=tellus_gfz_datafile,
                      tellus_scalefile=tellus_scalefile,
                      grace_dates_file=grace_dates_file, gracetype=gracetype,
                      max_distance=max_distance, datasource=datasource, greensfunc=greensfunc, scaling=scaling);
    print(MyParams);
    print("\n");
    return MyParams;
