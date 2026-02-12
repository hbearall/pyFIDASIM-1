# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 11:24:43 2021

@author: bgeiger3
"""
import numpy as np
import copy

def calc_arot_brot(direction):
    assert(direction.size == 3)
    y = direction[2]
    x = np.sqrt(np.sum(direction[:]**2))
    b = np.arctan2(y, x)
    Arot = np.array([[np.cos(b), 0., np.sin(b)], [
                    0., 1., 0.], [-np.sin(b), 0., np.cos(b)]])
    y = direction[1]
    x = direction[0]
    a = np.arctan2(y, x) - np.pi
    Brot = np.array([[np.cos(a), -np.sin(a), 0.],
                     [np.sin(a), np.cos(a), 0.], [0., 0., 1.]])
    return Arot, Brot


def nbi_geometry():
    """Fills a blank dictionary with nbi geometry related parameters.
    
    Returns
    ----------
    nbi_geometry : dictionary, nbi geometry parameters
    """
    
    default_geometry = {
        'ion_source_size': np.array([22.8, 50.6]), 
        'focal_length': np.array([650., 700.]), 
        'divergence': np.array([0.8, 0.8]) / 180 * np.pi, 
        'aperture_1_rectangular' : True, 
        'aperture_1_size': np.array([100., 100.]), 
        'aperture_1_distance': 650., 
        'aperture_1_offset': np.array([0., 0.]), 
        'aperture_2_rectangular' : True, 
        'aperture_2_size': np.array([100., 100.]), 
        'aperture_2_distance': 650.,   # cm
        'aperture_2_offset': np.array([0., 0.])
    }
    
    Q1_geometry = copy.deepcopy(default_geometry)
    Q1_geometry['ID'] = 'Q1'
    Q1_geometry['source_position'] = np.array(
        [607.5594, 1171.7889, 29.5000])  # cm
    Q1_geometry['direction'] = np.array([-0.365400, -0.926946, -0.085174])
    Arot, Brot = calc_arot_brot(Q1_geometry['direction'])
    Q1_geometry['Arot'] = Arot
    Q1_geometry['Brot'] = Brot

    
    Q2_geometry = copy.deepcopy(default_geometry)
    Q2_geometry['ID'] = 'Q2'
    Q2_geometry['source_position'] = np.array(
        [692.2962, 1131.0989, 29.5000])  # cm
    Q2_geometry['direction'] = np.array([-0.494953, -0.864735, -0.085174])
    Arot, Brot = calc_arot_brot(Q2_geometry['direction'])
    Q2_geometry['Arot'] = Arot
    Q2_geometry['Brot'] = Brot
    
    
    Q3_geometry = copy.deepcopy(default_geometry)
    Q3_geometry['ID'] = 'Q3'
    Q3_geometry['source_position'] = np.array(
        [692.2962, 1131.0989, -90.5000])  # cm
    Q3_geometry['direction'] = np.array([-0.494953, -0.864735, 0.085174])
    Arot, Brot = calc_arot_brot(Q3_geometry['direction'])
    Q3_geometry['Arot'] = Arot
    Q3_geometry['Brot'] = Brot
    
    
    Q4_geometry = copy.deepcopy(default_geometry)
    Q4_geometry['ID'] = 'Q4'
    Q4_geometry['source_position'] = np.array(
        [607.5594, 1171.7889, -90.5000])  # cm
    Q4_geometry['direction'] = np.array([-0.365400, -0.926946, 0.085174])
    Arot, Brot = calc_arot_brot(Q4_geometry['direction'])
    Q4_geometry['Arot'] = Arot
    Q4_geometry['Brot'] = Brot
    
    
    Q5_geometry = copy.deepcopy(default_geometry)
    Q5_geometry['ID'] = 'Q5'
    Q5_geometry['source_position'] = np.array(
        [197.2344, 1305.1116, -29.5000])  # cm
    Q5_geometry['direction'] = np.array([-0.249230, -0.964692, 0.085174])
    Arot, Brot = calc_arot_brot(Q5_geometry['direction'])
    Q5_geometry['Arot'] = Arot
    Q5_geometry['Brot'] = Brot
    
    
    
    Q6_geometry = copy.deepcopy(default_geometry)
    Q6_geometry['ID'] = 'Q6'
    Q6_geometry['source_position'] = np.array(
        [104.7639, 1321.9997, -29.5000])  # cm
    Q6_geometry['direction'] = np.array([-0.107854, -0.990511, 0.085174])
    Arot, Brot = calc_arot_brot(Q6_geometry['direction'])
    Q6_geometry['Arot'] = Arot
    Q6_geometry['Brot'] = Brot
    
    
    Q7_geometry = copy.deepcopy(default_geometry)
    Q7_geometry['ID'] = 'Q7'
    Q7_geometry['source_position'] = np.array([104.27, 1317.11, 90.67])
    #from David: np.array([104.76, 1321.998, 90.5])#from pyfidasim gitlab: np.array([104.27, 1317.11, 90.67])  # cm
    Q7_geometry['direction'] = np.array([-0.1085, -0.99, -.085])
    #from David: np.array([-0.10785, -0.9905, -0.0862])#from pyfidasim gitlab: np.array([-0.1085, -0.99, -.085])
    Arot, Brot = calc_arot_brot(Q7_geometry['direction'])
    Q7_geometry['Arot'] = Arot
    Q7_geometry['Brot'] = Brot
    Q7_geometry['aperture_1_distance']= 626.9   # cm
    Q7_geometry['aperture_1_offset'] = np.array([1.7,-4.8])
    Q7_geometry['aperture_1_size']   = np.array([34.1, 67.0])
    
    Q7_geometry['aperture_2_distance']= 708.4   # cm
    Q7_geometry['aperture_2_offset'] = np.array([1.9,3.2])
    Q7_geometry['aperture_2_size']   = np.array([34.1, 57.4])
    
    
    Q8_geometry = copy.deepcopy(default_geometry)
    Q8_geometry['ID'] = 'Q8'
    Q8_geometry['source_position'] = np.array([195.96, 1300.37, 90.67])  # cm
    Q8_geometry['direction'] = np.array([-0.2486, -0.9648, -.085])
    Q8_geometry['aperture_1_distance']= 626.9   # cm
    Q8_geometry['aperture_1_offset'] = np.array([5.6,-4.8])
    Q8_geometry['aperture_1_size']   = np.array([34.1, 67.0])
    
    Q8_geometry['aperture_2_distance']= 708.4   # cm
    Q8_geometry['aperture_2_offset'] = np.array([-5.9,3.2])
    Q8_geometry['aperture_2_size']   = np.array([34.1, 57.4])
    
    Arot, Brot = calc_arot_brot(Q8_geometry['direction'])
    Q8_geometry['Arot'] = Arot
    Q8_geometry['Brot'] = Brot
    
    nbi_geometry = {}
    nbi_geometry['sources'] = ['Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8']
    nbi_geometry['Q1'] = Q1_geometry
    nbi_geometry['Q2'] = Q2_geometry
    nbi_geometry['Q3'] = Q3_geometry
    nbi_geometry['Q4'] = Q4_geometry
    nbi_geometry['Q5'] = Q5_geometry
    nbi_geometry['Q6'] = Q6_geometry
    nbi_geometry['Q7'] = Q7_geometry
    nbi_geometry['Q8'] = Q8_geometry
    
    return nbi_geometry