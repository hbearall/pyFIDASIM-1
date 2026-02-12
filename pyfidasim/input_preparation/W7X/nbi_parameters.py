# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 11:24:43 2021

@author: bgeiger3
"""
#----------------------------------------------------------------------------------
#- routines to read NBI parameter information for a given shot and temporal range -
#----------------------------------------------------------------------------------

# General Note:
#   While the routine itself works from any directory, saving the data only works
#   when called from the 'examples/W7X' folder where the typical start file is located
#   For bug report and necessary fixes/alterations: thir@ipp.mpg.de

# Input:
#   shot_number: shot number in W7-X standard format YYYYMMDD.XXX
#   t_start:     starting time of NBI operation window
#   t_stop:      end time of NBI operation window
#   fractions:   variable to, if demanded, pass fractions which are then used instead of the ones written to the file
#
#   General structure is suited to assesses all NBI operation in a single shot
#   by passing t_start and t_stop as arrays
#
# example input nbi_parameters('20180920.042', t_start = 6.5, t_stop = 6.52, debug = True)

# OUTPUT:
#   beam_information:  dictionary suited as a pyFIDASIM input for the beam parameters - source resolved
#                      also save to the 'examples/W7X/Data/' folder

# Import general packages
import os
import numpy as np
import matplotlib.pyplot as plt

# import the routines for loading & saving hdf5 files
import pyfidasim
# from pyfidasim.hdf5 import save_dict, load_dict
from pyfidasim.toolbox import save_dict, load_dict


###############################################################################
# routine to read the NBI parameters from the W7-X data base & a lookup table #
###############################################################################
def nbi_parameters(shot_number = '20180920.042', t_start = 6.5, t_stop = 6.52, fractions = None, debug = False,default=False):
    defaultFractions = [0.3, 0.5, 0.2]
    if default:
        import copy
        nbi_params = {}
        nbi_params['sources'] = ['Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8']
        
        params = {}
        params['voltage'] = 54.6  # keV
        params['power'] = 1.5  # MW
        params['fraction'] = [0.3, 0.5, 0.2]
        params['type'] = 'H'
    
        for source in nbi_params['sources']:
            nbi_params[source] = copy.deepcopy(params)
        return nbi_params
    if fractions == 'default':
        fractions = defaultFractions
    
    # import 'custom' W7-X routines
    from w7xdia.nbi     import get_nspec_currentFractions
    # THE FOLLOWING IMPORTS NEED TO BE EXPANDED WITH ADDITIONAL SOURCES BECOMING AVAILABLE IN THE FUTURE
    from w7xdia.nbi import get_source_7_power,   get_source_8_power
    from w7xdia.nbi import get_source_7_voltage, get_source_8_voltage

    # put together the variable shot information - checks whether passed times are an array or single floats
    # float case    
    try:
        len(t_start)
        # turn single floats into arrays
        shot_information = [shot_number, t_start, t_stop]
    # array case
    except:
        # just put the variables together without altering them
        shot_information = [shot_number, [t_start], [t_stop]]
    
    # set a time out value for the database 
    timeout = 10

    # set up a dictionary with the routines - NEEDS TO BE EXPANDED WITH NEW SOURCES
    routines = {}
    # sub-levels for the voltage, energy and beam fractions
    routines['voltage']       = {}
    routines['voltage']['Q7'] = get_source_7_voltage
    routines['voltage']['Q8'] = get_source_8_voltage
    routines['power']         = {}
    routines['power'] ['Q7']  = get_source_7_power
    routines['power'] ['Q8']  = get_source_8_power
    routines['fractions']     = get_nspec_currentFractions
    # set up a dictionary to contain general source informations as well as the information which is later saved time resolved
    source_information                 = {}
    source_information['sources']      = ['Q7','Q8']#['Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8']
    source_information['source_index'] = [6,7]      #[0,1,2,3,4,5,6,7]

    # array of all the active sources and the respective indices
    sources      = []
    source_index = []
    # threshold power value to determine whether a source is No documentation available on or not
    power_threshold = 0.1   # MW
    volt_threshold  = 40    # keV

    #########################
    # loop over all sources #
    #########################
    for idx, sour in enumerate(source_information['sources']):
        
        # set up a sub-dictionary for the source
        source_information[sour] = {}

        #############################################################
        # load source power from data base outside of temporal loop #
        #############################################################
        # as one import contains all the power and voltage data for a given shot,
        # load them before the different time stamps are assessed individually
        source_information[sour]['time_power'], source_information[sour]['power'] = routines['power'][sour](shot_information[0], timeout = timeout)

        # assign a source the 'turned on'-status if it's output power surpasses the given threshold
        if max(source_information[sour]['power']) > power_threshold:
            sources.append(sour)
            source_index.append(source_information['source_index'][idx])

        #####################################################################################
        # calculate the average output power per beam by averaging the passed time interval #
        #####################################################################################
        # only do that if the source was identified as relavant beforehand however
        if sour in sources:
            # now also loop over all passed time ranges
            for PIT_idx, PIT in enumerate(shot_information[1]):
                # add a sub level for the point in time to the dictionary
                handle    = 't = ' + str(PIT)
                source_information[sour][handle] = {}
                # rename variables to have shorter names
                time      = source_information[sour]['time_power']
                power     = source_information[sour]['power']
                # reduce the power array based on the time slot
                power_red = power[np.where(time <= shot_information[2][PIT_idx])]
                time_red  = time [np.where(time <= shot_information[2][PIT_idx])]
                power_red = power_red[np.where(time_red >= PIT)]
                time_red  = time_red [np.where(time_red >= PIT)]
                # do the averaging math
                power_avg = sum(power_red) / len(power_red)
                # save the result to the dictionary
                source_information[sour][handle]['power'] = power_avg

            # if wanted, plot the source time trend
            if debug:
                plt.plot(time, power, label = 'Source ' + str(source_information['sources'][idx]) + '[MW]')
                plt.xlim(left=0)
                plt.xlabel('t in s')
                plt.ylabel('NBI power in MW')
                plt.legend()
                plt.title('%s - Power trace' % (shot_information[0]))
                plt.show()

        #####################################################################
        # read the acceleration voltage, again outside of the temporal loop #
        #####################################################################
        source_information[sour]['time_voltage'], source_information[sour]['voltage'] = routines['voltage'][sour](shot_information[0], timeout = timeout)

        # if the source was previously identified to be relevant, assess its voltage
        if sour in sources:
            # now also loop over all the passed time points
            for PIT_idx, PIT in enumerate(shot_information[1]):
                # put together the handle for the dictionary
                handle   = 't = ' + str(PIT)
                # rename variables to have shorter names
                time     = source_information[sour]['time_voltage']
                voltage  = source_information[sour]['voltage']
                # reduce the arrays for the relevant range
                volt_red = voltage[(time <= shot_information[2][PIT_idx]) & (time >= shot_information[1][PIT_idx])]
                time_red = time   [(time <= shot_information[2][PIT_idx]) & (time >= shot_information[1][PIT_idx])]
                # reduce the data points to only those above the threshold as we cannot afford to have 
                # an averaged voltage with zeros as we need the applied voltage
                volt_red = volt_red[(volt_red > volt_threshold)]
                # do the averaging math
                volt_avg = np.average(volt_red)
                # save the result to the dictionary
                source_information[sour][handle]['voltage'] = volt_avg

            # plot the data if demanded
            if debug:
                plt.plot(time, voltage*1.e-3, label = 'Source ' + str(source_information['sources'][idx]) + '[kV]')
                plt.xlim(left=0)
                plt.xlabel('t in s')
                plt.ylabel('NBI voltage in kV')
                plt.legend()
                plt.title('%s - Voltage trace' % (shot_information[0]))
                plt.show()

        ##################################################
        # read the beam fractions from the look up table #
        ##################################################

        # assess all the different sources
        if sour in sources:
            # read in the beam fractions for the given source only when None are passed
            if fractions == None:
                try:
                    t, sig, unc = get_nspec_currentFractions(shot_number, sour)
                    indrange = np.where((t >= shot_information[1]) & (t <= shot_information[2]))[0]
                    beam_fractions = []
                    for ii in range(sig.shape[0]):
                        beam_fractions = np.append(beam_fractions, np.mean(sig[ii, indrange]))
                    print(beam_fractions)
                    add_shot_beam_fractions(shot_number, beam_fractions, sour, overwrite = False)
                except ValueError:
                    beam_fractions = read_beam_fractions(shot_number, sour)
            else:
                beam_fractions = fractions
            # now also loop over all the NBI active points in time
            for PIT_idx, PIT in enumerate(shot_information[1]):
                # put together the handle for the dictionary
                handle   = 't = ' + str(PIT)
                # currently implemented via look up - needs to be changed in the long term
                source_information[sour][handle]['fraction'] = beam_fractions
                #get_nspec_currentFractions(shot_information[0], source_information['source_index'][idx])

        ###############################################################################
        # eventually write some general information which is the same for all sources #
        ###############################################################################
        source_information[sour]['type'] = 'H'    # obviously hardcoded atm, may change later
        source_information['source_arr'] = source_index

    # print the as used identified sources
    print('Sources %s found to be active in discharge' % sources)
        
    # get the path to the storage folder
    path = os.path.dirname(pyfidasim.__file__)
    file_path = os.path.abspath(os.path.join(path,os.pardir))
        
    # save the dictionary for all the points in time resolved to individual hdf5 files
    for PIT_idx, PIT in enumerate(shot_information[1]):
        # set up a dictionary for the data related to a single point in time
        beam_information = {}
        # put together the handle for the dictionary
        handle   = 't = ' + str(PIT)
        # fill the dictionary with some general information
        beam_information['sources']    = source_information['sources']
        beam_information['source_arr'] = source_information['source_arr']

        # loop over the sources
        for sour in source_information['sources']:
            # set up a sub-dictionary
            beam_information[sour] = {}
            # if the source was identified to be active earlier, save the data
            if sour in sources:
                # fill the dictionary with source specific information
                beam_information[sour]['voltage']    = source_information[sour][handle]['voltage']
                beam_information[sour]['power']      = source_information[sour][handle]['power']
                beam_information[sour]['fraction']   = source_information[sour][handle]['fraction']
                beam_information[sour]['type']       = source_information[sour]['type']
            # if the source was identified to be inactive, fill the dictionary with some generic values
            else:
                beam_information[sour]['voltage']    = 0
                beam_information[sour]['power']      = 0
                beam_information[sour]['fraction']   = [1.0, 0.0, 0.0] 
                beam_information[sour]['type']       = 'H'
                
            # get rid of nans in case of alternating sources in same shot
            voltage = np.array([beam_information[sour]['voltage'] ])
            power = np.array([beam_information[sour]['power'] ])
            
            if np.isnan(voltage) or voltage <=0.0: beam_information[sour]['voltage'] = 0.0
            if np.isnan(power) or power <= 0.0: beam_information[sour]['power'] = 0.0
                
        # put together the file name and the path to the data folder
        # I tried making the path nicer but the other package caused unnecessary complications
        file_path += '/examples/W7X/Data/'
        file_path = os.path.abspath(file_path)
        file_name = str(shot_information[0]) + '_' + str(int(round(PIT*1000))) + '_NBI_parameters.hdf5'
        # and save the dictionary as an hdf5 file - in try beacuse does work from any other directory
        # try:
        #     save_dict(beam_information, file_path + os.sep + file_name)
        # except OSError:
        #     print('Failed at writing file due to assumed path structure')
            
    return beam_information

#############################################################
# routine to read the beam fractions from a hdf5 dictionary #
#############################################################
def read_beam_fractions(shot_number, source):
    # shot number specifies the day the data was taken at
    # as parameters seems to change on a day to day base no time information within the shot asked for
    
    # put together the path to the .hdf5 file
    path_name = os.path.dirname(pyfidasim.__file__)
    path_name = os.path.dirname(path_name)
    path_name += '/examples/W7X/Data/'
    path_name = os.path.abspath(path_name)
    # path_name = os.path.abspath(os.path.join(path_name, os.pardir))
    
    # load the dictionary
    data = load_dict(path_name + os.sep + 'W7X_beam_fractions_lookup.hdf5')
    
    # check whether data for the requested shot number is available
    try:
        beam_fractions  = data[shot_number][source]
        return beam_fractions        
    except:
        print('\nNo beam fractions for requested shot %s found in look up dictionary' % (shot_number))
        print('To add data, call add_shot_beam_fractions - path to data is:')
        print('http://archive-webapi.ipp-hgw.mpg.de/ArchiveDB/raw/W7XAnalysis/CDX_NI_Spectroscopy/USB_HR4000_Q8_NeutralParticleFractions_DATASTREAM/V2/0/E1')
        print()
        raise ValueError('No beam fractions stored for shot %s' % (shot_number))
    
################################################################################
# routine to add an entry to the hdf5 dictionary containing the beam fractions #
################################################################################
# example input:
# add_shot_beam_fractions('20180920.042', [0.31, 0.54, 0.15], 'Q8')

def add_shot_beam_fractions(shot_number, beam_fractions, source, overwrite = False):
    
    # check for the inputs to have proper format - crude but should find obvious mistakes
    if len(beam_fractions) != 3:
        raise ValueError('Length of beam fraction array not applicable')
    day, program = shot_number.split('.')
    if len(day) != 8 or len(program) != 3:
        raise ValueError('Passed shot number not in correct format - needs to be yyyymmdd.xxx')
    try:
        int(day)
        int(program)
    except:
        raise ValueError('Passed shot number has non-number character')
    
    # put together the path to the .hdf5 file
    path_name = os.path.dirname(pyfidasim.__file__)
    path_name = os.path.abspath(os.path.join(path_name, os.pardir))
    path_name += '/examples/W7X/Data/'
    file_name = 'W7X_beam_fractions_lookup.hdf5'
    
    # open the dictionary
    data = load_dict(path_name + file_name)
    
    # check whether the requested shot already exists
    if shot_number in data:
        # check whether data for the given source is already present
        if source in data[shot_number].keys():
            if not overwrite:
                print('\nShot %s with fractions %s already present in dictionary' % (shot_number, beam_fractions))
                print('To overwrite data, call again with overwrite = True\n')
            elif overwrite:
                print('\nOverwriting fractions for shot %s to %s\n' % (shot_number, beam_fractions))
                data[shot_number][source] = beam_fractions
        else:
            data[shot_number][source] = beam_fractions
            print('\nAdded %s as fractions for shot %s\n' % (beam_fractions, shot_number))
    else:
        data[shot_number] = {}
        data[shot_number][source] = beam_fractions
        print('\nAdded %s as fractions for shot %s\n' % (beam_fractions, shot_number))
        
    # save the dict with the additional entry
    save_dict(data, path_name + file_name)
    
##############
# test input #
##############
if __name__ == '__main__':
    # one source example
    beam_information = nbi_parameters('20180920.042', t_start = 6.5, t_stop = 6.52, fractions=[0.3,0.5,0.2], debug = True)
    # two source example
    beam_information = nbi_parameters('20180918.046', t_start = 2.7, t_stop = 3.9, debug = True)
    pass