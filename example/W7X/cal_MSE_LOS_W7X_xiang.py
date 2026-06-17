### Calculation for MSE-LOS (2x8 fibre array) from AEA21 port, red shifts
### Code for pyFIDASIM run on W7X. This pyFIDASIM is a new version developed by Aiden.
 
import numpy as np

from pyfidasim.input_preparation.W7X.los_geometry import los_geometry
from pyfidasim.input_preparation.W7X.equilibrium import get_reference_equilibrium#, equilibrium
from pyfidasim.input_preparation.W7X import nbi#, plasma_profiles
import matplotlib.pyplot as plt
import functionsBES as bes

# import matplotlib as matplotlib
# matplotlib.rcParams.update({'font.size': 22.5})
# plt.close('all')
colors = plt.get_cmap('Dark2')

expID = '20250515.028'#'20180919.039'
tcal = 3.7#[3.55, 3.6]
LOS_head = 'MSE_LOS'#'los_BES'#'AEA21_A:'

wavelen_range = [652., 662.] # nm
drz_phi_ran = [2., 0.45*np.pi, 0.55*np.pi]
nmarker = int(10.e3)

''' define the sim_settings '''
fidasim = {}
fidasim["nmarker"] = nmarker
fidasim['calc_halo'] = True
fidasim['separate_dcx'] = True # If True, separates Direct Charge Exchange (DCX) from other halo contributions in output.
fidasim['verbose'] = True
fidasim['calc_photon_origin'] = False
fidasim['calc_photon_origin_type'] = False
fidasim['calc_PSF'] = False
# fidasim['image_pos'] = bes.LOS_MSE_op24()['los_on_beam']
# fidasim['image_vec'] = bes.LOS_MSE_op24()['los_vec']
fidasim['calc_spectra'] = True
fidasim['calc_density'] = True # If True, stores beam/halo neutral density on the standard machine grid (R, Z, Phi)
fidasim['calc_rzp_dens'] = True
fidasim['calc_uvw'] = True
fidasim['calc_extended_emission'] = False
fidasim['respawun_if_aperture_is_hit'] = True
fidasim['seed'] = nmarker//10 # Random number generator seed. Set to -1 for random runs, or a fixed integer for reproducibility
fidasim["batch_marker"] = fidasim["nmarker"]

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
# fidasim['nbi_sources'] = 'Q7'

from pyfidasim.input_prep import input_prep
sim_settings,spec,tables,fields,profiles,nbi,grid3d,ncdf,PSF,fbm = input_prep(fidasim)
# Placeholders
# ncdf = {'active': False, 'lambda0': 656.1, 'trans': 1, 'l_to_dwp': 1, 'spectrum_extended': False, 'cdfvars': np.zeros((1,1,1))}
# PSF = {'image_pos': 1, 'image_vec': 1, 'image_blur': 0, 'los_image_arr': 0, 'n_rand': 0, 'f_lens': 1}
# fbm = {'afbm': 0, 'fbm': np.zeros((1,1,1,1)), 'denf': np.zeros((1,1,1)), 'btipsign': -1, 'emin':0,'eran':0,'nenergy':0,'energy':np.zeros(1),'dE':0,'pmin':0,'pran':0,'npitch':0,'pitch':np.zeros(1),'dP':0}

if 'denimp' not in profiles.keys(): # impurity density
    zimp = tables['zimps'][0]
    profiles['denp'] = np.zeros((len(profiles['s'])))
    profiles['denimp'] = np.zeros((1,len(profiles['s'])))
    profiles['denimp'][0,:]=profiles['dene'] * (profiles['zeff'] - 1) / (zimp**2 - zimp)
    profiles['denp']=profiles['dene'] - profiles['denimp'][0,:]*zimp
# -------------------------------
# -- Plot LOS and NBI -----------
# -------------------------------

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
grid3d,spec = calc_attenuation(sim_settings,profiles,nbi,spec,fields,grid3d,tables,ncdf,PSF,fbm)
t2 = time.time()

print("Time taken: %.2f minutes."%((t2-t1)/60.))

# -----------------------------
# Plot Beam Emission Spectra
# -----------------------------
if 'intens' in spec:
    fig_spec, ax_slider = plot_spectra_interactive(spec, labels=['full', 'half', 'third', 'halo'])
    fig_spec.show()
# ax_slider.set_visible(False) # temporally hide the sliding bar to save the figure

# -----------------------------
# Plot Radial Density Profiles
# -----------------------------
fig_dens, ax_den = plot_midplane_heatmap(grid3d, fields, energy_indices=[0, 1, 2])
for i, k in enumerate(spec['losname']):
    xyzk = spec['los_pos'][i,] + spec['full_los_length'][i]*spec['los_vec'][i,]
    ax_den.plot([spec['los_pos'][i,0], xyzk[0]], [spec['los_pos'][i,1], xyzk[1]], 
                marker = '', lw = 1., color = colors(i%4), label = k)


from pyfidasim.los import nbi_intersection
rzp_los, distances = nbi_intersection(spec, nbi)
xlos = rzp_los[:,:,0] * np.cos(rzp_los[:,:,2])
ylos = rzp_los[:,:,0] * np.sin(rzp_los[:,:,2])
ax_den.plot(xlos, ylos, marker = 'o', mfc = None, color = 'k', lw = 0.0, markersize = 3,
            alpha = 0.5)
# -----------------------------
## plot s-along LOS combined with the NBI density distribution along the LOS
from pyfidasim.los import calc_s_along_los
### plot the LoS on the beam radiance's contour
# for ilos in range(spec['nlos']):
#     dist_arr,xyz_arr,s_arr = calc_s_along_los(spec['los_pos'][ilos,:],spec['los_vec'][ilos,:],fields,dl=0.1)
#     ax_den.plot(xyz_arr[0,], xyz_arr[1,], lw = 0.4, color = colors(ilos%len(colors.colors)))
# fig_dens.show()

nr_plot=160
npos_max=len(spec['s_per_grid_intersection'][0,:])
density_along_los = np.zeros((spec['nlos'],npos_max))
distance_along_los_arr = np.empty_like(density_along_los)

ra_grid=np.linspace(0,1.2,num=nr_plot)
resolution=np.zeros((spec['nlos'],nr_plot))
# ax1 = plt.figure().add_subplot(projection = '3d')
for ilos in range(spec['nlos']):
    for ipos in range(npos_max):
        ## define NBI density along LOS
        index=spec['los_grid_intersection_indices'][ilos,ipos,:]
        if not np.all(index == 0):
            density_along_los[ilos,ipos]=np.nansum(grid3d['density'][:,:,index[0],index[1],index[2]],axis=(0,1))
        
        ## fill the resolution plot
        ra=np.sqrt(spec['s_per_grid_intersection'][ilos, ipos])
        ii=np.argmin(np.abs(ra_grid-ra))
        resolution[ilos,ii]+=density_along_los[ilos,ipos]
    ## Normalize the density along LOS
    # density_along_los[ilos,:]/=np.max(density_along_los[ilos,:]) 
    ## Normalize the resolution data
    resolution[ilos,:]/=np.max(resolution[ilos,:])
    
    dl = spec['dl_per_grid_intersection'][ilos, :]
    s = spec['s_per_grid_intersection'][ilos, :]
    distance_along_los_arr[ilos,] = spec['full_los_length'][ilos]-np.sum(dl)+np.cumsum(dl)-0.5*dl
    
plt.figure()
for ilos in range(spec['nlos']):
    idx = density_along_los[ilos,] > 0
    dist_arr,xyz_arr,s_arr = calc_s_along_los(spec['los_pos'][ilos,:],spec['los_vec'][ilos,:],fields,dl=0.1)
    ra_val = np.min(s_arr)
    
    plt.fill_between(distance_along_los_arr[ilos, idx], np.zeros(np.sum(idx)) + ra_val,
                     ra_val + density_along_los[ilos,idx]*0.1,
                     color=colors(ilos%len(colors.colors)), alpha=0.3)

plt.xlabel('distance along LOS [cm]')
plt.ylabel('r/a')
# plt.ylim([0,1])
plt.tight_layout()
    


