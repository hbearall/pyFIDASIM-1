### Code for pyFIDASIM run on W7X. This pyFIDASIM is a new version developed by Aiden.
 
import numpy as np

from pyfidasim.input_preparation.W7X.los_geometry import los_geometry
from pyfidasim.input_preparation.W7X.equilibrium import get_reference_equilibrium#, equilibrium
from pyfidasim.input_preparation.W7X import nbi#, plasma_profiles
# import matplotlib.pyplot as plt
# import matplotlib as matplotlib
# matplotlib.rcParams.update({'font.size': 22.5})
# plt.close('all')

shotnum = '20180919.039'
twin_profile = [3.55, 3.6]
LOS_head = 'AEA21_A:0'

wavelen_range = [652., 662.] # nm
drz_phi_ran = [2., 0.45*np.pi, 0.55*np.pi]
nmarker = int(10.e3)

''' define the sim_settings '''
fidasim = {}
fidasim["nmarker"] = nmarker
fidasim['calc_halo'] = True
fidasim['verbose'] = True
fidasim['calc_photon_origin'] = True
fidasim['calc_photon_origin_type'] = True
fidasim['calc_PSF'] = False
fidasim['calc_spectra'] = True
fidasim['calc_density'] = True # If True, stores beam/halo neutral density on the standard machine grid (R, Z, Phi)
fidasim['calc_rzp_dens'] = True
fidasim['calc_uvw'] = True
fidasim['calc_extended_emission'] = False
fidasim['respawun_if_aperture_is_hit'] = True
fidasim['seed'] = 12345 # Random number generator seed. Set to -1 for random runs, or a fixed integer for reproducibility
fidasim["batch_marker"] = fidasim["nmarker"]
fidasim['separate_dcx'] = False # If True, separates Direct Charge Exchange (DCX) from other halo contributions in output.

''' setup for W7X run'''
fidasim['progID'] = shotnum
fidasim['t_start'] = twin_profile[0]
fidasim['t_stop'] = twin_profile[1]
fidasim['machine'] = 'W7X'

''' setup the spectrometer '''
fidasim['los_head'] = LOS_head
fidasim['los_file'] = '//share.ipp-hgw.mpg.de/documents/xiha/Documents/Python_W7X/pyFIDASIM/example/Data/op12b_ils_geometry.txt'
fidasim['los_default'] = True
fidasim['los_new'] = False

spec = los_geometry(head = fidasim['los_head'], file = fidasim['los_file'], default = True, new = False)
fidasim['los_pos'] = spec['los_pos']
fidasim['los_vec'] = spec['los_vec']
fidasim['los_name'] = spec['losname']
fidasim['only_pi'] = False
fidasim['dlam'] = 1.e-2 # nm
fidasim['lambda_min'] = wavelen_range[0]
fidasim['lambda_max'] = wavelen_range[1]
# fidasim['output_individual_stark_lines'] = True
del spec

'''setup the equilibrium'''
vmecID, b0_scaling = get_reference_equilibrium(fidasim['progID'])
fidasim['vmecID'] = vmecID
fidasim['b0_factor'] = np.round(b0_scaling, 4)
# fidasim['drz'] = drz_phi_ran[0]

''' setup the grid '''
# grid3d = define_grid3d(sgrid, drz=1., r_ran = [520., 620.], z_ran = [5., 55.], phi_ran=[0.45*np.pi,0.55*np.pi])
fidasim['grid_drz'] = drz_phi_ran[0] # precision of the grid in cm
fidasim['r_ran'] = [520., 620.] #major radius extent of the simulation domain [cm]
fidasim['z_ran'] = [5., 55.] # vertical extent of the simulation domain [cm]
fidasim['phi_ran'] = drz_phi_ran[1:] # toroidal angle extent [radians]
fidasim['u_range'] = [] #[700,830] # extent of the beam-aligned grid along the beam propagation axis [cm]
# fidasim['u_in_vessel'] = True # added by Xiang, re-define the u_range to include the beam in the vessel only. The trimesh package is needed.
fidasim['v_width'] = 60.0 # total width of the beam-aligned grid perpendicular to beam (horizontal) [cm]
fidasim['w_width'] = 30.0 # total width of the beam-aligned grid perpendicular to beam (vertical) [cm]
fidasim['du'] = 1. # step length along beam direction (u) for beam grid [cm]
fidasim['dvw'] = 1.e-1 # step length perpendicular to beam (v, w) for beam grid [cm]

fidasim['impurities'] = ['Carbon']
# fidasim['path_to_tables'] = ''
# fidasim['load_raw_data'] = True

fidasim['lambda_min'] = wavelen_range[0]
fidasim['lambda_max'] = wavelen_range[1]
fidasim['dlam'] = 1.e-2 # nm

''' define the nbi parameters '''
fidasim['shot_num'] = shotnum
fidasim['nbi_default'] = False
fidasim['nbi_debug'] = False
fidasim['nbi_mass'] = 1. # 1 Hydrogen
fidasim['ion_mass'] = 1. # 1 for protium
fidasim['nbi_sources'] = 'Q7'

from pyfidasim.input_prep import input_prep
sim_settings,spec,tables,fields,profiles,nbi,grid3d,ncdf,PSF,fbm = input_prep(fidasim)
# Placeholders
ncdf = {'active': False, 'lambda0': 656.1, 'trans': 1, 'l_to_dwp': 1, 'spectrum_extended': False, 'cdfvars': np.zeros((1,1,1))}
PSF = {'image_pos': 1, 'image_vec': 1, 'image_blur': 0, 'los_image_arr': 0, 'n_rand': 0, 'f_lens': 1}
fbm = {'afbm': 0, 'fbm': np.zeros((1,1,1,1)), 'denf': np.zeros((1,1,1)), 'btipsign': -1, 'emin':0,'eran':0,'nenergy':0,'energy':np.zeros(1),'dE':0,'pmin':0,'pran':0,'npitch':0,'pitch':np.zeros(1),'dP':0}

if 'denimp' not in profiles.keys(): # impurity density
    zimp = tables['zimps'][0]
    profiles['denp'] = np.zeros((len(profiles['s'])))
    profiles['denimp'] = np.zeros((1,len(profiles['s'])))
    profiles['denimp'][0,:]=profiles['dene'] * (profiles['zeff'] - 1) / (zimp**2 - zimp)
    profiles['denp']=profiles['dene'] - profiles['denimp'][0,:]*zimp
# profiles['denimp'] *= 1.e1
# -------------------------------
# -- Plot LOS and NBI -----------
# -------------------------------
'''
fields['Er'], 
'''
from pyfidasim.plotting_routines import plot_geometry_3d, plot_spectra_interactive, plot_midplane_heatmap, plot_profiles, plot_magnetic_equilibrium
fig_w7x_geo = plot_geometry_3d(fields, spec, nbi=nbi, grid3d=grid3d, plot_crossed_cells = False)
fig_w7x_geo.show()
fig_w7x_mag = plot_magnetic_equilibrium(fields)
fig_profiles = plot_profiles(profiles)
# ----------------------------
# -- run pyFIDASIM -----------
# ----------------------------
import time
from pyfidasim.main import calc_attenuation

t1 = time.time()
grid3d,spec=calc_attenuation(sim_settings,profiles,nbi,spec,fields,grid3d,tables,ncdf,PSF,fbm)
t2 = time.time()

print("Time taken: ",(t2-t1))

# -----------------------------
# Plot Beam Emission Spectra
# -----------------------------
if 'intens' in spec:
    fig_spec = plot_spectra_interactive(spec, labels=['full', 'half', 'third', 'halo'])
    fig_spec.show()

# -----------------------------
# Plot Radial Density Profiles
# -----------------------------
fig_dens = plot_midplane_heatmap(grid3d, fields)
fig_dens.show()
# -----------------------------
# Plot Beam-Path Density Profiles
# -----------------------------
# from pyfidasim.plotting_routines import plot_u_density
# plot_u_density(grid3d)


