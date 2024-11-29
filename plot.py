#[ ........................................................... ]#
#[
#[ Maxwell Rosen
#[ plot.py
#[ 2021-04-20
#  For plotting comparison of uniform and nonuniform grid runs.
#[
#[ ........................................................... ]#

import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from matplotlib.colors import LogNorm
import multiprocessing
from scipy.integrate import cumulative_trapezoid as cumtrapz
# import imageio.v2 as imageio

# dataDir = '/home/mr1884/scratch/Link to scratch_traverse/gkylmax/traverse-wham1x-compare_unif_vs_nonunif/outputs/'
# dataDir = './data-hires-lorad/'
dataDir = './'
mc2pFolder = ''
mc2pFilename = 'gk_wham_1x2v_p1_uniform'
mc2pUniformFolder = ''
mc2pUniformFilename = 'gk_step_1x2v_p1_cons_nonuniform'
numericFolder = ''
numericFilename = 'gk_step_1x2v_p1_cons_nonuniform'
frame_max_plus1 = 1
time_per_frame = 1e-6

plot_potential_trace = 0
plot_bimax_moms = 1
plot_subtracted_moms = 0
plot_integrate_positivity = 0

# frame_arr = np.arange(0,11)
# frame_arr = np.array([1:4])
# save_figure_as_file= 1     #[ If True, save figure to file. If False, display figure on screen.

# plot_moments       = 1  # Plot density, potential, upar, tperp, tpar.
# plot_distvpar      = 0  # plot distribution function in vpar.
# plot_distmu        = 0  # plot distribution function in mu.
# plot_distf_at_z    = 0
# z_loctions = np.array([0, 0.3, 0.98, 2.4])

# print(frame_arr)
# def process_frame(frameNum):
outDir = './'

figureFileFormat = '.png'    #[ Can be .png, .pdf, .ps, .eps, .svg.

#   #[ ............... End of user inputs (MAYBE) ..................... ]#

polyOrder = 1
basisType = 'gkhyb'

eps0, mu0 = 8.8541878176204e-12, 1.2566370614359e-06
eV        = 1.602176487e-19
qe, qi    = -1.602176487e-19, 1.602176487e-19
me, mp    = 9.10938215e-31, 1.672621637e-27

mi        = 2.014*mp                         #[ Deuterium ion mass.
Te0       = 940*eV
n0        = 3.e19
B_p       = 0.53
beta      = 0.4                              #[ Ratio of plasma to magnetic pressure.
tau       = (B_p**2)*beta/(2*mu0*n0*Te0)-1    #[ Ti/Te ratio.
Ti0       = tau*Te0
n_pol     = 3e19

timestep = 8.42244e-12

#[ Thermal speeds.
vti = np.sqrt(Ti0/mi)
vte = np.sqrt(Te0/me)
mui0 = 0.5*mi*(vti**2)/B_p
c_s = np.sqrt(Te0/mi)

#[ Gyrofrequencies and gyroradii.
omega_ci = eV*B_p/mi
rho_s    = c_s/omega_ci

#[ Some RGB colors. These are MATLAB-like.
defaultBlue    = [0, 0.4470, 0.7410]
defaultOrange  = [0.8500, 0.3250, 0.0980]
defaultGreen   = [0.4660, 0.6740, 0.1880]
defaultPurple  = [0.4940, 0.1840, 0.5560]
defaultRed     = [0.6350, 0.0780, 0.1840]
defaultSkyBlue = [0.3010, 0.7450, 0.9330]
grey           = [0.5, 0.5, 0.5]
#[ Colors in a single array.
defaultColors = [defaultBlue,defaultOrange,defaultGreen,defaultPurple,defaultRed,defaultSkyBlue,grey,'black']

#[ LineStyles in a single array.
lineStyles = ['-','--',':','-.','None','None','None','None']
markers    = ['None','None','None','None','o','d','s','+']

#[ Some fontsizes used in plots.
xyLabelFontSize       = 17
titleFontSize         = 17
colorBarLabelFontSize = 17
tickFontSize          = 14
legendFontSize        = 14
textFontSize          = 16

#.Set the font size of the ticks to a given size.
def setTickFontSize(axIn,fontSizeIn):
  axIn.tick_params(labelsize = fontSizeIn)

#.Plot vertical lines at +/- given x location.
def plot_verticalLinesPM(xIn, axIn):
  ymin, ymax = axIn.get_ylim()
  eps = 0.5*ymax,
  axIn.plot([xIn, xIn], [ymin-eps, ymax+eps], linestyle=":", color='grey')
  axIn.plot([-xIn, -xIn], [ymin-eps, ymax+eps], linestyle=":", color='grey')
  axIn.set_ylim(ymin, ymax)

#   #................................................................................#

if plot_potential_trace:
  print("Plotting potential trace")
  filename_bmag = str(dataDir+mc2pFilename+'-bmag.gkyl')
  pgData_bmag = pg.GData(filename_bmag)
  pgInterp_bmag = pg.GInterpModal(pgData_bmag, polyOrder, 'ms')
  x_bmag, dataOut_bmag = pgInterp_bmag.interpolate()

  bmag_shape = dataOut_bmag.shape
  midpoint = int(bmag_shape[0]/2)
  upperhalf = dataOut_bmag[midpoint:]
  peak = np.argmax(upperhalf)
  peak_idx = midpoint+peak

  def loadphi(frame_number, filename):
    filename_phi = str(dataDir+filename+'-field_'+str(frame_number)+'.gkyl')
    pgData_phi = pg.GData(filename_phi)
    pgInterp_phi = pg.GInterpModal(pgData_phi, polyOrder, 'ms')
    x_phi, dataOut_phi = pgInterp_phi.interpolate()
    return dataOut_phi
  
  def get_temp(frame_number, filename):
    # filename_elc = str(dataDir+filename+'-elc_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
    # pgData_elc = pg.GData(filename_elc)
    # pgInterp_elc = pg.GInterpModal(pgData_elc, polyOrder, 'ms')
    # coords, Tpar_elc = pgInterp_elc.interpolate(2)
    # coords, Tperp_elc = pgInterp_elc.interpolate(3)
    # Temp = (Tpar_elc[midpoint,0] + 2*Tperp_elc[midpoint,0])/3 * me / eV
    return Te0

  potential = np.zeros(frame_max_plus1)
  potential_mod = np.zeros(frame_max_plus1)
  Temp = np.zeros(frame_max_plus1)
  Temp_mod = np.zeros(frame_max_plus1)  
  for i in range(frame_max_plus1):
    dataOut_phi = loadphi(i, mc2pFilename)
    Temp[i] = get_temp(i, mc2pFilename)
    midphi = dataOut_phi[midpoint]
    phi_peak = dataOut_phi[peak_idx]
    potential[i] = (midphi[0] - phi_peak[0]) / Temp[i]

    dataOut_phi = loadphi(i, numericFilename)
    Temp_mod[i] = get_temp(i, numericFilename)
    midphi = dataOut_phi[midpoint]
    phi_peak = dataOut_phi[peak_idx]
    potential_mod[i] = (midphi[0] - phi_peak[0]) / Temp_mod[i]

  Temp = Temp*eV
  Temp_mod = Temp_mod*eV

  plt.plot(np.arange(frame_max_plus1)*1e-6, potential, label = 'Standard collisions')
  plt.plot(np.arange(frame_max_plus1)*1e-6, potential_mod, linestyle='--', label = 'Modified collisions')
  plt.xlabel('Time, seconds')
  plt.ylabel('Potential difference, e phi / T_e(psi_{min},z=0)')
  plt.title('Potential difference between midplane and peak magnetic field')
  plt.legend()
  plt.savefig(outDir+'potential_trace'+figureFileFormat)
  plt.close()

  plt.plot(np.arange(frame_max_plus1)*1e-6, Temp/eV, label = 'Standard collisions')
  plt.plot(np.arange(frame_max_plus1)*1e-6, Temp_mod/eV, linestyle='--', label = 'Modified collisions')
  plt.xlabel('Time, seconds')
  plt.ylabel('Temperature, T_e(psi_{min},z=0)')
  plt.title('Temperature at midplane')
  plt.legend()
  plt.savefig(outDir+'temperature_trace'+figureFileFormat)
  plt.close()

if plot_bimax_moms:
  def make_moms(frame_number):
    print("Getting moments for frame ", frame_number)
    # filename_elc = str(dataDir+unifFile+'-elc_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
    # pgData_elc = pg.GData(filename_elc)
    # pgInterp_elc = pg.GInterpModal(pgData_elc, polyOrder, 'ms')
    # coords, n_elc = pgInterp_elc.interpolate(0)
    # coords, u_elc = pgInterp_elc.interpolate(1)
    # coords, Tpar_elc = pgInterp_elc.interpolate(2)
    # coords, Tperp_elc = pgInterp_elc.interpolate(3)

    filename_ion_mc2p = str(mc2pFolder+mc2pFilename+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
    pgData_ion_mc2p = pg.GData(filename_ion_mc2p)
    pgInterp_ion_mc2p = pg.GInterpModal(pgData_ion_mc2p, polyOrder, 'ms')
    coords, n_ion_mc2p = pgInterp_ion_mc2p.interpolate(0)
    coords, u_ion_mc2p = pgInterp_ion_mc2p.interpolate(1)
    coords, Tpar_ion_mc2p = pgInterp_ion_mc2p.interpolate(2)
    coords, Tperp_ion_mc2p = pgInterp_ion_mc2p.interpolate(3)

    filename_field_mc2p = str(mc2pFolder+mc2pFilename+'-field_'+str(frame_number)+'.gkyl')
    pgData_field_mc2p = pg.GData(filename_field_mc2p)
    pgInterp_field_mc2p = pg.GInterpModal(pgData_field_mc2p, polyOrder, 'ms')
    coords, phi_mc2p = pgInterp_field_mc2p.interpolate()

    filename_ion_mc2p_uniform = str(mc2pUniformFolder+''+mc2pUniformFilename+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
    pgData_ion_mc2p_uniform = pg.GData(filename_ion_mc2p_uniform)
    pgInterp_ion_mc2p_uniform = pg.GInterpModal(pgData_ion_mc2p_uniform, polyOrder, 'ms')
    coords, n_ion_mc2p_uniform = pgInterp_ion_mc2p_uniform.interpolate(0)
    coords, u_ion_mc2p_uniform = pgInterp_ion_mc2p_uniform.interpolate(1)
    coords, Tpar_ion_mc2p_uniform = pgInterp_ion_mc2p_uniform.interpolate(2)
    coords, Tperp_ion_mc2p_uniform = pgInterp_ion_mc2p_uniform.interpolate(3)

    filename_field_mc2p_uniform = str(mc2pUniformFolder+''+mc2pUniformFilename+'-field_'+str(frame_number)+'.gkyl')
    pgData_field_mc2p_uniform = pg.GData(filename_field_mc2p_uniform)
    pgInterp_field_mc2p_uniform = pg.GInterpModal(pgData_field_mc2p_uniform, polyOrder, 'ms')
    coords, phi_mc2p_uniform = pgInterp_field_mc2p_uniform.interpolate()

    # filename_elc_mod = str(dataDir+modifiedFile+'-elc_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
    # print("Reading file ", filename_elc_mod)
    # pgData_elc_mod = pg.GData(filename_elc_mod)
    # pgInterp_elc_mod = pg.GInterpModal(pgData_elc_mod, polyOrder, 'ms')
    # coords, n_elc_mod = pgInterp_elc_mod.interpolate(0)
    # coords, u_elc_mod = pgInterp_elc_mod.interpolate(1)
    # coords, Tpar_elc_mod = pgInterp_elc_mod.interpolate(2)
    # coords, Tperp_elc_mod = pgInterp_elc_mod.interpolate(3)
    
    filename_ion_numeric = str(numericFolder+''+numericFilename+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
    pgData_ion_numeric = pg.GData(filename_ion_numeric)
    pgInterp_ion_numeric = pg.GInterpModal(pgData_ion_numeric, polyOrder, 'ms')
    coords, n_ion_numeric = pgInterp_ion_numeric.interpolate(0)
    coords, u_ion_numeric = pgInterp_ion_numeric.interpolate(1)
    coords, Tpar_ion_numeric = pgInterp_ion_numeric.interpolate(2)
    coords, Tperp_ion_numeric = pgInterp_ion_numeric.interpolate(3)
    
    filename_field_numeric = str(numericFolder+''+numericFilename+'-field_'+str(frame_number)+'.gkyl')
    pgData_field_numeric = pg.GData(filename_field_numeric)
    pgInterp_field_numeric = pg.GInterpModal(pgData_field_numeric, polyOrder, 'ms')
    coords, phi_numeric = pgInterp_field_numeric.interpolate()

    data = pg.GData(str(mc2pFolder+''+mc2pFilename+"-nodes.gkyl"))
    vals = data.get_values()
    nodes_R_mc2p = vals[:,0]
    nodes_Z_mc2p = vals[:,1]
    nodes_phi_mc2p = vals[:,2]

    data_mc2p_uniform = pg.GData(str(mc2pUniformFolder+''+mc2pUniformFilename+"-nodes.gkyl"))
    vals_mc2p_uniform = data_mc2p_uniform.get_values()
    nodes_R_mc2p_uniform = vals_mc2p_uniform[:,0]
    nodes_Z_mc2p_uniform = vals_mc2p_uniform[:,1]
    nodes_phi_mc2p_uniform = vals_mc2p_uniform[:,2]

    data_numeric = pg.GData(str(numericFolder+''+numericFilename+"-nodes.gkyl"))
    vals_numeric = data_numeric.get_values()
    nodes_R_numeric = vals_numeric[:,0]
    nodes_Z_numeric = vals_numeric[:,1]
    nodes_phi_numeric = vals_numeric[:,2]

    shape_R = np.shape(nodes_R_mc2p)
    midplane_R_min = nodes_R_mc2p[shape_R[0]//2]
    midplane_R_max = nodes_R_mc2p[shape_R[0]//2]
    throat_R_min = nodes_R_mc2p[shape_R[0]//4]
    throat_R_max = nodes_R_mc2p[shape_R[0]//4]

    def expand_1D_array(original_array):
      new_length = 2 * len(original_array) - 1
      new_array = np.zeros(new_length)
      new_array[0] = original_array[0]
      for i in range(1, len(original_array)):
          new_array[2*i - 1] = (original_array[i - 1] + original_array[i]) / 2
          new_array[2*i] = original_array[i]
      return new_array
    
    def expand_2D_array(original_array):
      original_shape = np.shape(original_array)
      new_shape = (2*original_shape[0]-1, 2*original_shape[1]-1)
      new_array = np.zeros(new_shape)
      for i in range(1, original_shape[0]-1):
        for j in range(1, original_shape[1]-1):
          new_array[2*i, 2*j] = original_array[i, j]
          new_array[2*i, 2*j+1] = (original_array[i, j] + original_array[i, j+1]) / 2
          new_array[2*i, 2*j-1] = (original_array[i, j] + original_array[i, j-1]) / 2
          new_array[2*i+1, 2*j] = (original_array[i, j] + original_array[i+1, j]) / 2
          new_array[2*i-1, 2*j] = (original_array[i, j] + original_array[i-1, j]) / 2
          new_array[2*i+1, 2*j+1] = (original_array[i, j] + original_array[i+1, j+1]) / 2
          new_array[2*i-1, 2*j-1] = (original_array[i, j] + original_array[i-1, j-1]) / 2
          new_array[2*i+1, 2*j-1] = (original_array[i, j] + original_array[i+1, j-1]) / 2
          new_array[2*i-1, 2*j+1] = (original_array[i, j] + original_array[i-1, j+1]) / 2
      new_array[:,0] = expand_1D_array(original_array[:,0])
      new_array[:,-1] = expand_1D_array(original_array[:,-1])
      new_array[0,:] = expand_1D_array(original_array[0,:])
      new_array[-1,:] = expand_1D_array(original_array[-1,:])
      return new_array
    
    def expand_3D_array(original_array):
      original_shape = np.shape(original_array)
      new_shape = (2*original_shape[0]-1, 2*original_shape[1]-1, 2*original_shape[2]-1)
      new_array = np.zeros(new_shape)
      for i in range(1, original_shape[0]-1):
        for j in range(1, original_shape[1]-1):
          for k in range(1, original_shape[2]-1):
            new_array[2*i, 2*j, 2*k] = original_array[i, j, k]

            new_array[2*i, 2*j, 2*k+1] = (original_array[i, j, k] + original_array[i, j, k+1]) / 2
            new_array[2*i, 2*j, 2*k-1] = (original_array[i, j, k] + original_array[i, j, k-1]) / 2
            new_array[2*i, 2*j+1, 2*k] = (original_array[i, j, k] + original_array[i, j+1, k]) / 2
            new_array[2*i, 2*j-1, 2*k] = (original_array[i, j, k] + original_array[i, j-1, k]) / 2
            new_array[2*i+1, 2*j, 2*k] = (original_array[i, j, k] + original_array[i+1, j, k]) / 2
            new_array[2*i-1, 2*j, 2*k] = (original_array[i, j, k] + original_array[i-1, j, k]) / 2

            new_array[2*i+1, 2*j+1, 2*k] = (original_array[i, j, k] + original_array[i+1, j+1, k]) / 2
            new_array[2*i-1, 2*j-1, 2*k] = (original_array[i, j, k] + original_array[i-1, j-1, k]) / 2
            new_array[2*i+1, 2*j-1, 2*k] = (original_array[i, j, k] + original_array[i+1, j-1, k]) / 2
            new_array[2*i-1, 2*j+1, 2*k] = (original_array[i, j, k] + original_array[i-1, j+1, k]) / 2

            new_array[2*i+1, 2*j, 2*k+1] = (original_array[i, j, k] + original_array[i+1, j, k+1]) / 2
            new_array[2*i-1, 2*j, 2*k-1] = (original_array[i, j, k] + original_array[i-1, j, k-1]) / 2
            new_array[2*i+1, 2*j, 2*k-1] = (original_array[i, j, k] + original_array[i+1, j, k-1]) / 2
            new_array[2*i-1, 2*j, 2*k+1] = (original_array[i, j, k] + original_array[i-1, j, k+1]) / 2

            new_array[2*i, 2*j+1, 2*k+1] = (original_array[i, j, k] + original_array[i, j+1, k+1]) / 2
            new_array[2*i, 2*j-1, 2*k-1] = (original_array[i, j, k] + original_array[i, j-1, k-1]) / 2
            new_array[2*i, 2*j+1, 2*k-1] = (original_array[i, j, k] + original_array[i, j+1, k-1]) / 2
            new_array[2*i, 2*j-1, 2*k+1] = (original_array[i, j, k] + original_array[i, j-1, k+1]) / 2
      new_array[:,0,:] = expand_2D_array(original_array[:,0,:])
      new_array[:,-1,:] = expand_2D_array(original_array[:,-1,:])
      new_array[:,:,0] = expand_2D_array(original_array[:,:,0])
      new_array[:,:,-1] = expand_2D_array(original_array[:,:,-1])
      new_array[0,:,:] = expand_2D_array(original_array[0,:,:])
      new_array[-1,:,:] = expand_2D_array(original_array[-1,:,:])
      return new_array

    
    nodes_Z_mc2p = expand_1D_array(nodes_Z_mc2p)
    nodes_Z_mc2p_uniform = expand_1D_array(nodes_Z_mc2p_uniform)
    nodes_Z_numeric = expand_1D_array(nodes_Z_numeric)
    nodes_Z_mc2p = nodes_Z_mc2p[1:]
    nodes_Z_mc2p_uniform = nodes_Z_mc2p_uniform[1:]
    nodes_Z_numeric = nodes_Z_numeric[1:]
    # nodes_R = expand_1D_array(nodes_R)


    # n_elc = n_elc[:,0]
    # u_elc = u_elc[:,0]
    # Tpar_elc = Tpar_elc[:,0] * me / eV
    # Tperp_elc = Tperp_elc[:,0] * me / eV
    # T_elc = (Tpar_elc + 2*Tperp_elc)/3
    n_ion_mc2p = n_ion_mc2p[:,0]
    u_ion_mc2p = u_ion_mc2p[:,0]
    Tpar_ion_mc2p = Tpar_ion_mc2p[:,0] * mi / eV
    Tperp_ion_mc2p = Tperp_ion_mc2p[:,0] * mi / eV
    T_ion_mc2p = (Tpar_ion_mc2p + 2*Tperp_ion_mc2p)/3
    phi_mc2p = phi_mc2p[:,0]
    # midplane_Te = T_elc[T_elc.shape[0]//2]
    ephioTe_mc2p =  qi * phi_mc2p / Te0

    n_ion_mc2p_uniform = n_ion_mc2p_uniform[:,0]
    u_ion_mc2p_uniform = u_ion_mc2p_uniform[:,0]
    Tpar_ion_mc2p_uniform = Tpar_ion_mc2p_uniform[:,0] * mi / eV
    Tperp_ion_mc2p_uniform = Tperp_ion_mc2p_uniform[:,0] * mi / eV
    T_ion_mc2p_uniform = (Tpar_ion_mc2p_uniform + 2*Tperp_ion_mc2p_uniform)/3
    phi_mc2p_uniform = phi_mc2p_uniform[:,0]
    # midplane_Te = T_elc[T_elc.shape[0]//2]
    ephioTe_mc2p_uniform =  qi * phi_mc2p_uniform / Te0

    # n_elc_mod = n_elc_mod[:,0]
    # u_elc_mod = u_elc_mod[:,0]
    # Tpar_elc_mod = Tpar_elc_mod[:,0] * me / eV
    # Tperp_elc_mod = Tperp_elc_mod[:,0] * me / eV
    # T_elc_mod = (Tpar_elc_mod + 2*Tperp_elc_mod)/3
    n_ion_numeric = n_ion_numeric[:,0]
    u_ion_numeric = u_ion_numeric[:,0]
    Tpar_ion_numeric = Tpar_ion_numeric[:,0] * mi / eV
    Tperp_ion_numeric = Tperp_ion_numeric[:,0] * mi / eV
    T_ion_numeric = (Tpar_ion_numeric + 2*Tperp_ion_numeric)/3
    phi_numeric = phi_numeric[:,0]
    # midplane_Te_mod = T_elc_mod[T_elc_mod.shape[0]//2]
    ephioTe_numeric =  qi * phi_numeric / Te0

    # # Compute polarization density for ions
    # # Read in the magnetic field
    # filename_bmag = str(dataDir+unifFile+'-bmag.gkyl')
    # pgData_bmag = pg.GData(filename_bmag)
    # pgInterp_bmag = pg.GInterpModal(pgData_bmag, polyOrder, 'ms')
    # coords, bmag = pgInterp_bmag.interpolate()
    # bmag_shape = bmag.shape

    # # Checked in gkyl_gk_geometry_bmag_mid that this is how it's done
    # epsilon_i = mi * n_pol / bmag[bmag_shape[0]//2,bmag_shape[1]//2,0]**2

    # filename_jacobgeo = str(dataDir+unifFile+'-jacobgeo.gkyl')
    # pgData_jacobgeo = pg.GData(filename_jacobgeo)
    # pgInterp_jacobgeo = pg.GInterpModal(pgData_jacobgeo, polyOrder, 'ms')
    # coords, jacobgeo = pgInterp_jacobgeo.interpolate()
    # jacobgeo = jacobgeo[:,0]

    # filename_gxx = str(dataDir+unifFile+'-gxxj.gkyl')
    # pgData_gxx = pg.GData(filename_gxx)
    # pgInterp_gxx = pg.GInterpModal(pgData_gxx, polyOrder, 'ms')
    # coords, gxxj = pgInterp_gxx.interpolate()
    # gxxj = gxxj[:,0]

    # D =  gxxj * epsilon_i
    # dpsi = coords[1][0] - coords[0][0]
    # ni_pol = np.zeros(D.shape)
    # for i in range (D.shape[0]-1):
    #   ni_pol[i,:] = -1/jacobgeo[i,:] / dpsi * (((D[i+1,:]+D[i,:])/2) * ((phi[i+1,:]-phi[i,:])/dpsi) - \
    #                                            ((D[i-1,:]+D[i,:])/2) * ((phi[i,:]-phi[i-1,:])/dpsi))
    # ni_pol[0,:] = 0.0

    # make an array grid that is the size of coords

    nonunif_mapc2p_filename = str(mc2pFolder+''+mc2pFilename+'-mapc2p.gkyl')
    pgData_mc2p_mapc2p = pg.GData(nonunif_mapc2p_filename)
    pgInterp_mc2p_mapc2p = pg.GInterpModal(pgData_mc2p_mapc2p, polyOrder, 'ms')
    x_nonunif_mapc2p, dataOut_mapc2p_mc2p = pgInterp_mc2p_mapc2p.interpolate(0)

    unif_mapc2p_filename = str(mc2pUniformFolder+''+mc2pUniformFilename+'-mapc2p.gkyl')
    pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
    pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
    x_unif_mapc2p, dataOut_mapc2p_mc2p_uniform = pgInterp_unif_mapc2p.interpolate(2)

    numeric_mapc2p_filename = str(numericFolder+''+numericFilename+'-mapc2p.gkyl')
    pgData_numeric_mapc2p = pg.GData(numeric_mapc2p_filename)
    pgInterp_numeric_mapc2p = pg.GInterpModal(pgData_numeric_mapc2p, polyOrder, 'ms')
    x_unif_mapc2p, dataOut_numeric_mapc2p = pgInterp_numeric_mapc2p.interpolate(2)


    X_mc2p = nodes_Z_mc2p
    X_mc2p_uniform = nodes_Z_mc2p_uniform
    X_numeric = nodes_Z_numeric

    # X = nodes_Z[:,:]
    # Y = nodes_R[:,:]
    
    fig, ax = plt.subplots(3, 3, figsize=(12,12))
    fig.suptitle(str(frame_number*time_per_frame)+' seconds', fontsize=20)

    def plot_moment_data(data_mc2p, data_mc2p_uniform, data_numeric, ax, fig, title, locx, locy):
      ax[locx,locy].plot(X_mc2p, data_mc2p, label=mc2pFilename, color='blue')
      ax[locx,locy].plot(X_mc2p_uniform, data_mc2p_uniform, label=mc2pUniformFilename, linestyle='-.', color='orange')
      ax[locx,locy].plot(X_numeric, data_numeric, label=numericFilename, linestyle='--', color='green')
      ax[locx,locy].set_xlabel('Z cylindrical axis, m')
      ax[locx,locy].set_ylabel(title)
      ax[locx,locy].set_title(title, fontsize=16)

    # plot_moment_data(n_elc, n_elc_mod, ax, fig, 'n_e, m^{-3}', 0, 0)
    # plot_moment_data(u_elc, u_elc_mod, ax, fig, 'U_{e,||}, m/s', 0, 2)
    # plot_moment_data(Tpar_elc, Tpar_elc_mod, ax, fig, 'T_{e,||}, eV', 1, 0)
    # plot_moment_data(Tperp_elc, Tperp_elc_mod, ax, fig, 'T_{e,perp}, eV', 1, 1)
    # plot_moment_data(T_elc, T_elc_mod, ax, fig, 'T_e, eV', 1, 2)

    plot_moment_data(n_ion_mc2p, n_ion_mc2p_uniform, n_ion_numeric, ax, fig, 'n_i, m^{-3}', 0, 0)
    plot_moment_data(u_ion_mc2p, u_ion_mc2p_uniform,u_ion_numeric, ax, fig, 'U_{i,||}, m/s', 0, 2)
    plot_moment_data(Tpar_ion_mc2p, Tpar_ion_mc2p_uniform, Tpar_ion_numeric, ax, fig, 'T_{i,||}, eV', 1, 0)
    plot_moment_data(Tperp_ion_mc2p, Tperp_ion_mc2p_uniform, Tperp_ion_numeric, ax, fig, 'T_{i,perp}, eV', 1, 1)
    plot_moment_data(T_ion_mc2p, T_ion_mc2p_uniform, T_ion_numeric, ax, fig, 'T_i, eV', 1, 2)

    ax[0,2].legend()

    # Plot electron density on a log scale
    # ax[0,1].plot(X,n_elc, label='Standard collisions')
    # ax[0,1].plot(X,n_elc_mod, label='Modified collisions', linestyle='--')
    # ax[0,1].set_yscale('log')
    # ax[0,1].set_xlabel('Z cylindrical axis, m')
    # ax[0,1].set_ylabel('n_e')
    # ax[0,1].set_title('n_e (log scale) m^{-3}', fontsize=16)

    # Plot the ion density on a log scale
    ax[0,1].plot(X_mc2p,n_ion_mc2p, label=mc2pFolder, color='blue')
    ax[0,1].plot(X_mc2p_uniform,n_ion_mc2p_uniform, label=mc2pUniformFolder, linestyle='-.', color='orange')
    ax[0,1].plot(X_numeric,n_ion_numeric, label=numericFolder, linestyle='--', color='green')
    ax[0,1].set_yscale('log')
    ax[0,1].set_xlabel('Z cylindrical axis, m')
    ax[0,1].set_ylabel('n_i')
    ax[0,1].set_title('n_i (log scale) m^{-3}', fontsize=16)

    plot_moment_data(phi_mc2p, phi_mc2p_uniform, phi_numeric, ax, fig, 'phi, V', 2, 0)
    plot_moment_data(ephioTe_mc2p, ephioTe_mc2p_uniform, ephioTe_numeric, ax, fig, 'e phi / T_e', 2, 1)
    ax[2,2].remove()

    plt.tight_layout()
    plt.savefig(outDir+'moments_'+str(frame_number)+figureFileFormat, dpi=600)
    plt.show()
    # plt.close()

    if plot_subtracted_moms:
      interp_x = np.linspace(X_mc2p[0], X_mc2p[-1], 1000)
      n_ion_interp_nunif = np.interp(interp_x, X_numeric, n_ion_numeric)
      n_ion_interp_unif = np.interp(interp_x, X_mc2p, n_ion_mc2p)
      n_ion_diff_rel = (n_ion_interp_unif - n_ion_interp_nunif) / n_ion_interp_unif
      n_ion_diff = n_ion_interp_unif - n_ion_interp_nunif

      u_ion_interp_nunif = np.interp(interp_x, X_numeric, u_ion_numeric)
      u_ion_interp_unif = np.interp(interp_x, X_mc2p, u_ion_mc2p)
      u_ion_diff_rel = (u_ion_interp_unif - u_ion_interp_nunif) / u_ion_interp_unif
      u_ion_diff = u_ion_interp_unif - u_ion_interp_nunif

      Tpar_ion_interp_nunif = np.interp(interp_x, X_numeric, Tpar_ion_numeric)
      Tpar_ion_interp_unif = np.interp(interp_x, X_mc2p, Tpar_ion_mc2p)
      Tpar_ion_diff_rel = (Tpar_ion_interp_unif - Tpar_ion_interp_nunif) / Tpar_ion_interp_unif
      Tpar_ion_diff = Tpar_ion_interp_unif - Tpar_ion_interp_nunif

      Tperp_ion_interp_nunif = np.interp(interp_x, X_numeric, Tperp_ion_numeric)
      Tperp_ion_interp_unif = np.interp(interp_x, X_mc2p, Tperp_ion_mc2p)
      Tperp_ion_diff_rel = (Tperp_ion_interp_unif - Tperp_ion_interp_nunif) / Tperp_ion_interp_unif
      Tperp_ion_diff = Tperp_ion_interp_unif - Tperp_ion_interp_nunif

      T_ion_interp_nunif = np.interp(interp_x, X_numeric, T_ion_numeric)
      T_ion_interp_unif = np.interp(interp_x, X_mc2p, T_ion_mc2p)
      T_ion_diff_rel = (T_ion_interp_unif - T_ion_interp_nunif) / T_ion_interp_unif
      T_ion_diff = T_ion_interp_unif - T_ion_interp_nunif

      phi_interp_nunif = np.interp(interp_x, X_numeric, phi_numeric)
      phi_interp_unif = np.interp(interp_x, X_mc2p, phi_mc2p)
      phi_diff_rel = (phi_interp_unif - phi_interp_nunif) / phi_interp_unif
      phi_diff = phi_interp_unif - phi_interp_nunif

      ephioTe_interp_nunif = np.interp(interp_x, X_numeric, ephioTe_numeric)
      ephioTe_interp_unif = np.interp(interp_x, X_mc2p, ephioTe_mc2p)
      ephioTe_diff_rel = (ephioTe_interp_unif - ephioTe_interp_nunif) / ephioTe_interp_unif
      ephioTe_diff = ephioTe_interp_unif - ephioTe_interp_nunif

      def plot_diff_data(data, ax, fig, title, locx, locy):
        ax[locx,locy].plot(interp_x, data)
        ax[locx,locy].set_xlabel('Z cylindrical axis, m')
        ax[locx,locy].set_ylabel(title)
        ax[locx,locy].set_title(title, fontsize=16)

      fig, ax = plt.subplots(4, 3, figsize=(12,16))
      fig.suptitle(str(frame_number*time_per_frame)+' seconds', fontsize=20)

      plot_diff_data(n_ion_diff, ax, fig, 'Delta n_i, m^{-3}', 0, 0)
      ax[0,0].set_ylim(-2e17, 2e17)
      plot_diff_data(u_ion_diff, ax, fig, 'Delta U_{i,||}, m/s', 0, 1)
      ax[0,1].set_ylim(-1e4, 1e4)
      plot_diff_data(Tpar_ion_diff, ax, fig, 'Delta T_{i,||}, eV', 0, 2)
      ax[0,2].set_ylim(-100, 100)
      plot_diff_data(Tperp_ion_diff, ax, fig, 'Delta T_{i,perp}, eV', 1, 0)
      ax[1,0].set_ylim(-100, 100)
      plot_diff_data(phi_diff, ax, fig, 'Delta phi, V', 1, 1)
      ax[1,1].set_ylim(-625, -550)
      plot_diff_data(ephioTe_diff, ax, fig, 'Delta e phi / T_e', 1, 2)
      ax[1,2].set_ylim(-0.70, -0.60)

      plot_diff_data(n_ion_diff_rel, ax, fig, 'Delta n_i / n_i', 2, 0)
      ax[2,0].set_ylim(-0.02, 0.01)
      plot_diff_data(u_ion_diff_rel, ax, fig, 'Delta U_{i,||} / U_{i,||}', 2, 1)
      ax[2,1].set_ylim(-0.1, 0.1)
      plot_diff_data(Tpar_ion_diff_rel, ax, fig, 'Delta T_{i,||} / T_{i,||}', 2, 2)
      ax[2,2].set_ylim(-0.1, 0.1)
      plot_diff_data(Tperp_ion_diff_rel, ax, fig, 'Delta T_{i,perp} / T_{i,perp}', 3, 0)
      ax[3,0].set_ylim(-0.05, 0.05)

      ax[3,1].remove()
      ax[3,2].remove()

      plt.tight_layout()
      plt.savefig(outDir+'moments_diff_'+str(frame_number)+figureFileFormat, dpi=600)
      plt.close()



  # Number of processes to run in parallel
  make_moms(0)




#   #....................................DEPRICATED CODE............................................#

#     densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_unif = pg.GData(densityFileName_unif)
#     pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
#     x_unif, dataOut_unif = pgInterp_unif.interpolate()


#     nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
#     pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
#     pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
#     x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)


#   if plot_moments:
#     dataOut_unif, dataOut_unif_mapc2p, dataOut_nonunif, dataOut_nonunif_mapc2p = load_mapped_data('-ion_M0_')
    
#     # Create a subfigure that is 2 by 3
#     fig, ax = plt.subplots(2, 3, figsize=(20,10))
#     # Plot the density
#     ax[0,0].plot(dataOut_unif_mapc2p[:,0], dataOut_unif[:,0],'r', label=unifFile)
#     ax[0,0].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif[:,0],'b--', label=nonunifFile)
#     ax[0,0].set_xlabel('Cylindrical length coordinate, Z (m)', fontsize=xyLabelFontSize)
#     ax[0,0].set_ylabel('n_i (m^{-3})', fontsize=xyLabelFontSize)
#     ax[0,0].legend(loc='upper left', fontsize=legendFontSize)
#     setTickFontSize(ax[0,0],tickFontSize)

#     # Plot phi
#     dataOut_unif, dataOut_unif_mapc2p, dataOut_nonunif, dataOut_nonunif_mapc2p = load_mapped_data('-field_')

#     dataOut_unif *= eV/Te0
#     dataOut_nonunif *= eV/Te0
#     # dataOut_reduced *= eV/Te0
#     # dataOut_coarse *= eV/Te0

#     ax[0,1].plot(dataOut_unif_mapc2p[:,0], dataOut_unif[:,0],'r')
#     ax[0,1].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif[:,0],'b--')
#     ax[0,1].set_xlabel('Cylindrical length coordinate, Z (m)', fontsize=xyLabelFontSize)
#     ax[0,1].set_ylabel('phi (V))', fontsize=xyLabelFontSize)
#     # ax[0,1].set_title('Frame '+str(frameNum))
#     ax[0,1].set_title("time = "+str(frameNum*1e-5)+" s")
#     setTickFontSize(ax[0,1],tickFontSize)

#     # Plot uPar
#     M0_unif, M0_map, M0_nonunif, M0_nonunif_map = load_mapped_data('-ion_M0_')
#     M1_unif, M1_map, M1_nonunif, M1_nonunif_map = load_mapped_data('-ion_M1_')

#     upar_unif = M1_unif[:,0]/M0_unif[:,0]
#     upar_nonunif = M1_nonunif[:,0]/M0_nonunif[:,0]

#     ax[1,0].plot(M0_map[:,0], upar_unif / c_s,'r')
#     ax[1,0].plot(M0_nonunif_map[:,0], upar_nonunif / c_s,'b--')
#     ax[1,0].set_xlabel('Cylindrical length coordinate, Z (m)', fontsize=xyLabelFontSize)
#     ax[1,0].set_ylabel('u_{parallel} / c_s (m/s)', fontsize=xyLabelFontSize)
#     setTickFontSize(ax[1,0],tickFontSize)
    
#     # Plot tPerp
#     # M0_unif, M0_map, M0_nonunif, M0_nonunif_map, M0_reduced, M0_reduced_map = load_mapped_data('-ion_M0_')
#     M2perp_unif, M2perp_map, M2perp_nonunif, M2perp_nonunif_map = load_mapped_data('-ion_M2perp_')

#     tPerp_unif = M2perp_unif[:,0]/M0_unif[:,0] * mi / eV
#     tPerp_nonunif = M2perp_nonunif[:,0]/M0_nonunif[:,0] * mi / eV

#     ax[1,1].plot(M0_map[:,0], tPerp_unif,'r')
#     ax[1,1].plot(M0_nonunif_map[:,0], tPerp_nonunif,'b--')
#     ax[1,1].set_xlabel('Cylindrical length coordinate, Z (m)', fontsize=xyLabelFontSize)
#     ax[1,1].set_ylabel('T_{perp} (eV)', fontsize=xyLabelFontSize)
#     setTickFontSize(ax[1,1],tickFontSize)

#     # Plot tPar
#     # M0_unif, M0_map, M0_nonunif, M0_nonunif_map, M0_reduced, M0_reduced_map = load_mapped_data('-ion_M0_')
#     # M1_unif, M1_map, M1_nonunif, M1_nonunif_map, M1_reduced, M1_reduced_map = load_mapped_data('-ion_M1_')
#     M2par_unif, M2par_map, M2par_nonunif, M2par_nonunif_map = load_mapped_data('-ion_M2par_')

#     tPar_unif = (M2par_unif[:,0] - M1_unif[:,0]**2/M0_unif[:,0]) * mi / eV / M0_unif[:,0]
#     tPar_nonunif = (M2par_nonunif[:,0] - M1_nonunif[:,0]**2/M0_nonunif[:,0]) * mi / eV / M0_nonunif[:,0]

#     ax[1,2].plot(M0_map[:,0], tPar_unif,'r')
#     ax[1,2].plot(M0_nonunif_map[:,0], tPar_nonunif,'b--')
#     ax[1,2].set_xlabel('Cylindrical length coordinate, Z (m)', fontsize=xyLabelFontSize)
#     ax[1,2].set_ylabel('T_{parallel} (eV)', fontsize=xyLabelFontSize)
#     setTickFontSize(ax[1,2],tickFontSize)

#     # Plot the grid and mapc2p
#     # ax[0,2].plot(dataOut_unif_mapc2p[:,0], dataOut_unif_mapc2p[:,0],'r', label='Uniform 280x96x192', markersize=0.5)
#     # ax[0,2].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif_mapc2p[:,0],'b.', label='Nonuniform 280x96x192', markersize=0.5)
#     # ax[0,2].plot(dataOut_reduced_mapc2p[:,0], dataOut_reduced_mapc2p[:,0],'g.', label='Nonuniform 140x96x192', markersize=0.5)
#     # ax[0,2].set_xlabel('Cylindrical length coordinate, Z (m)', fontsize=xyLabelFontSize)
#     # ax[0,2].set_ylabel('Mapped cylindrical length coordinate, Z (m)', fontsize=xyLabelFontSize)
#     # ax[0,2].legend(loc='upper left', fontsize=legendFontSize)
#     # setTickFontSize(ax[0,2],tickFontSize)
    
#     # print("mapc2p uniform grid: ", dataOut_unif_mapc2p[:,0])
#     # print("mapc2p nonuniform grid: ", dataOut_nonunif_mapc2p[:,0])
#     # print("last cell spacings uniform grid: ", dataOut_unif_mapc2p[-1,0] - dataOut_unif_mapc2p[-2,0])
#     # print("last cell spacings nonuniform grid: ", dataOut_nonunif_mapc2p[-1,0] - dataOut_nonunif_mapc2p[-2,0])
#     # print("ratio of last cell spacings: ", (dataOut_unif_mapc2p[-1,0] - dataOut_unif_mapc2p[-2,0])/(dataOut_nonunif_mapc2p[-1,0] - dataOut_nonunif_mapc2p[-2,0]))
    
#     figName = 'moments_'+str(frameNum)
#     if save_figure_as_file:
#       fig.savefig(outDir+figName+figureFileFormat)
#       plt.close()
#       print('Saved figure as '+figName+figureFileFormat)  
#     else:
#       plt.show()
  
#   if plot_distvpar:
#     # f_unif, f_map, f_nonunif, f_nonunif_map = load_mapped_data('-ion-')
#     dataName = '-ion_'
#     densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_unif = pg.GData(densityFileName_unif)
#     pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'gkhyb')
#     x_unif, dataOut_unif = pgInterp_unif.interpolate()
#     dataOut_unif = np.squeeze(dataOut_unif)

#     densityFileName_nonunif = str(dataDir+nonunifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_nonunif = pg.GData(densityFileName_nonunif)
#     pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'gkhyb')
#     x_nonunif, dataOut_nonunif = pgInterp_nonunif.interpolate()
#     dataOut_nonunif = np.squeeze(dataOut_nonunif)

#     nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
#     pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
#     pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
#     x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)
#     dataOut_nonunif_mapc2p = np.squeeze(dataOut_nonunif_mapc2p)

#     unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
#     pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
#     pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
#     x_unif_mapc2p, dataOut_unif_mapc2p = pgInterp_unif_mapc2p.interpolate(2)
#     dataOut_unif_mapc2p = np.squeeze(dataOut_unif_mapc2p)

#     # Convert from cell center to edges
#     zmin = x_unif_mapc2p[0][0]
#     zmax = x_unif_mapc2p[0][-1]
#     diffs  = dataOut_nonunif_mapc2p[0:-1] + np.diff(dataOut_nonunif_mapc2p)/2
#     edged_dataOut_nonunif_mapc2p = np.insert(diffs, 0, zmin)
#     edged_dataOut_nonunif_mapc2p = np.append(edged_dataOut_nonunif_mapc2p, zmax)

#     # Convert from cell center to edges
#     diffs  = dataOut_unif_mapc2p[0:-1] + np.diff(dataOut_unif_mapc2p)/2
#     edged_dataOut_unif_mapc2p = np.insert(diffs, 0, zmin)
#     edged_dataOut_unif_mapc2p = np.append(edged_dataOut_unif_mapc2p, zmax)

#     unif_Jgeo_filename = str(dataDir+unifFile+'-jacobgeo.gkyl')
#     pgData_unif_Jgeo = pg.GData(unif_Jgeo_filename)
#     pgInterp_unif_Jgeo = pg.GInterpModal(pgData_unif_Jgeo, polyOrder, 'ms')
#     x_unif_Jgeo, dataOut_unif_Jgeo = pgInterp_unif_Jgeo.interpolate()
#     dataOut_unif_Jgeo = np.squeeze(dataOut_unif_Jgeo)

#     nonunif_Jgeo_filename = str(dataDir+nonunifFile+'-jacobgeo.gkyl')
#     pgData_nonunif_Jgeo = pg.GData(nonunif_Jgeo_filename)
#     pgInterp_nonunif_Jgeo = pg.GInterpModal(pgData_nonunif_Jgeo, polyOrder, 'ms')
#     x_nonunif_Jgeo, dataOut_nonunif_Jgeo = pgInterp_nonunif_Jgeo.interpolate()
#     dataOut_nonunif_Jgeo = np.squeeze(dataOut_nonunif_Jgeo)

#     unif_distf_shape = dataOut_unif.shape
#     nonunif_distf_shape = dataOut_nonunif.shape

#     tile_unif_Jgeo = np.ones((unif_distf_shape[0], unif_distf_shape[1], unif_distf_shape[2]))
#     for i in range(unif_distf_shape[0]):
#         tile_unif_Jgeo[i,:,:] *= dataOut_unif_Jgeo[i]

#     tile_nonunif_Jgeo = np.ones((nonunif_distf_shape[0], nonunif_distf_shape[1], nonunif_distf_shape[2]))
#     for i in range(nonunif_distf_shape[0]):
#         tile_nonunif_Jgeo[i,:,:] *= dataOut_nonunif_Jgeo[i]

#     dataOut_unif = np.trapz(np.abs(dataOut_unif / tile_unif_Jgeo), axis=2)
#     dataOut_nonunif = np.trapz(np.abs(dataOut_nonunif / tile_nonunif_Jgeo), axis=2)

#     # Interpolate the non-uniform data onto a uniform grid
#     dataOut_unif_shape = dataOut_unif.shape
#     dataOut_nonunif_interp = np.zeros((dataOut_unif_shape[0], dataOut_unif_shape[1]))
#     for i in range(nonunif_distf_shape[1]):
#         dataOut_nonunif_interp[:,i] = np.interp(dataOut_unif_mapc2p, dataOut_nonunif_mapc2p, dataOut_nonunif[:,i])
#     data_difference = np.abs((dataOut_unif - dataOut_nonunif_interp))
#   #   # Need velocity space grids or to convert to edges of z for plotting

#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,12))

#     norm = LogNorm(vmin = 1e-14, vmax = 1e-4)  # Create a LogNorm instance

#     pcolormesh1 = ax1.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[1]*np.sqrt(2), dataOut_unif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot

#     # Label the axes
#     ax1.set_ylabel('vpar / vpar_max')
#     ax1.set_xlabel('Uniform grid Z, cylindrical coodinate (m)')
#     ax1.set_title('Frame '+str(frameNum))

#     # print(x_unif[1]/vti)
#     # print(x_nonunif[1]/vti)

#     # pcolormesh2 = ax2.pcolormesh(edged_dataOut_nonunif_mapc2p, x_nonunif[1]/vti, dataOut_nonunif.T, cmap='inferno', norm=norm)
#     pcolormesh2 = ax2.pcolormesh(edged_dataOut_nonunif_mapc2p, mapc2p_vel_vpar(x_nonunif[1]), dataOut_nonunif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
#     # Label the axes
#     ax2.set_ylabel('vpar / vpar_max')
#     ax2.set_xlabel('Nonuniform grid Z, cylindrical coodinate (m)') 
#   #   pcolormesh3 = ax3.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[1] / vti, data_difference.T, cmap='inferno', norm=norm)
#   #   fig.colorbar(pcolormesh3, ax=ax3)  # Add a colorbar to the plot
#   # # Label the axes
#   #   ax3.set_ylabel('vpar / vti')
#   #   ax3.set_xlabel('Absolute difference between distribution functions Z, cylindrical coordinate (m)')
    
#     figName = 'distf_vpar_'+str(frameNum)
#     if save_figure_as_file:
#       fig.savefig(outDir+figName+figureFileFormat)
#       plt.close()
#       print('Saved figure as '+figName+figureFileFormat)  
#     else:
#       plt.show()  

#   if plot_distmu:
#     # f_unif, f_map, f_nonunif, f_nonunif_map = load_mapped_data('-ion-')
#     dataName = '-ion_'
#     densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_unif = pg.GData(densityFileName_unif)
#     pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
#     x_unif, dataOut_unif = pgInterp_unif.interpolate()
#     dataOut_unif = np.squeeze(dataOut_unif)

#     densityFileName_nonunif = str(dataDir+nonunifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_nonunif = pg.GData(densityFileName_nonunif)
#     pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'ms')
#     x_nonunif, dataOut_nonunif = pgInterp_nonunif.interpolate()
#     dataOut_nonunif = np.squeeze(dataOut_nonunif)

#     nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
#     pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
#     pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
#     x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)
#     dataOut_nonunif_mapc2p = np.squeeze(dataOut_nonunif_mapc2p)

#     unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
#     pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
#     pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
#     x_unif_mapc2p, dataOut_unif_mapc2p = pgInterp_unif_mapc2p.interpolate(2)
#     dataOut_unif_mapc2p = np.squeeze(dataOut_unif_mapc2p)

#     # Convert from cell center to edges
#     zmin = x_unif_mapc2p[0][0]
#     zmax = x_unif_mapc2p[0][-1]
#     diffs  = dataOut_nonunif_mapc2p[0:-1] + np.diff(dataOut_nonunif_mapc2p)/2
#     edged_dataOut_nonunif_mapc2p = np.insert(diffs, 0, zmin)
#     edged_dataOut_nonunif_mapc2p = np.append(edged_dataOut_nonunif_mapc2p, zmax)

#     # Convert from cell center to edges
#     diffs  = dataOut_unif_mapc2p[0:-1] + np.diff(dataOut_unif_mapc2p)/2
#     edged_dataOut_unif_mapc2p = np.insert(diffs, 0, zmin)
#     edged_dataOut_unif_mapc2p = np.append(edged_dataOut_unif_mapc2p, zmax)

#     unif_Jgeo_filename = str(dataDir+unifFile+'-jacobgeo.gkyl')
#     pgData_unif_Jgeo = pg.GData(unif_Jgeo_filename)
#     pgInterp_unif_Jgeo = pg.GInterpModal(pgData_unif_Jgeo, polyOrder, 'ms')
#     x_unif_Jgeo, dataOut_unif_Jgeo = pgInterp_unif_Jgeo.interpolate()
#     dataOut_unif_Jgeo = np.squeeze(dataOut_unif_Jgeo)

#     nonunif_Jgeo_filename = str(dataDir+nonunifFile+'-jacobgeo.gkyl')
#     pgData_nonunif_Jgeo = pg.GData(nonunif_Jgeo_filename)
#     pgInterp_nonunif_Jgeo = pg.GInterpModal(pgData_nonunif_Jgeo, polyOrder, 'ms')
#     x_nonunif_Jgeo, dataOut_nonunif_Jgeo = pgInterp_nonunif_Jgeo.interpolate()
#     dataOut_nonunif_Jgeo = np.squeeze(dataOut_nonunif_Jgeo)

#     unif_distf_shape = dataOut_unif.shape
#     nonunif_distf_shape = dataOut_nonunif.shape

#     tile_unif_Jgeo = np.ones((unif_distf_shape[0], unif_distf_shape[1], unif_distf_shape[2]))
#     for i in range(unif_distf_shape[0]):
#         tile_unif_Jgeo[i,:,:] *= dataOut_unif_Jgeo[i]

#     tile_nonunif_Jgeo = np.ones((nonunif_distf_shape[0], nonunif_distf_shape[1], nonunif_distf_shape[2]))
#     for i in range(nonunif_distf_shape[0]):
#         tile_nonunif_Jgeo[i,:,:] *= dataOut_nonunif_Jgeo[i]

#     dataOut_unif = np.trapz(np.abs(dataOut_unif / tile_unif_Jgeo), axis=1)
#     dataOut_nonunif = np.trapz(np.abs(dataOut_nonunif / tile_nonunif_Jgeo), axis=1)

#     # Interpolate the non-uniform data onto a uniform grid
#     dataOut_unif_shape = dataOut_unif.shape
#     dataOut_nonunif_shape = dataOut_nonunif.shape
#     dataOut_nonunif_interp = np.zeros((dataOut_unif_shape[0], dataOut_unif_shape[1]))
#     for i in range(dataOut_nonunif_shape[1]):
#         dataOut_nonunif_interp[:,i] = np.interp(dataOut_unif_mapc2p, dataOut_nonunif_mapc2p, dataOut_nonunif[:,i])
#     data_difference = np.abs((dataOut_unif - dataOut_nonunif_interp))

#     # Need velocity space grids or to convert to edges of z for plotting
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,12))

#     norm = LogNorm(vmin = 1e-15, vmax = 1e-6)  # Create a LogNorm instance

#     pcolormesh1 = ax1.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[2], dataOut_unif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot
#     # Label the axes
#     ax1.set_ylabel('mu / mu_max')
#     ax1.set_xlabel('Uniform grid Z, cylindrical coodinate (m)')
#     ax1.set_title('Frame '+str(frameNum))

#     pcolormesh2 = ax2.pcolormesh(edged_dataOut_nonunif_mapc2p, mapc2p_vel_mu(x_nonunif[2]), dataOut_nonunif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
#     # Label the axes
#     ax2.set_ylabel('mu / mu_max')
#     ax2.set_xlabel('Nonuniform grid Z, cylindrical coodinate (m)')

#   #   pcolormesh3 = ax3.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[2] / vti, data_difference.T, cmap='inferno', norm=norm)
#   #   fig.colorbar(pcolormesh3, ax=ax3)  # Add a colorbar to the plot
#   # # Label the axes
#   #   ax3.set_ylabel('vpar / vti')
#   #   ax3.set_xlabel('Absolute difference between distribution functions Z, cylindrical coordinate (m)')
    
#     figName = 'distf_mu_'+str(frameNum)
#     if save_figure_as_file:
#       fig.savefig(outDir+figName+figureFileFormat)
#       plt.close()
#       print('Saved figure as '+figName+figureFileFormat)  
#     else:
#       plt.show() 
      
#   def plot_distf_at_z_eq(z0_coordinate):
#     dataName = '-ion_'
#     densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_unif = pg.GData(densityFileName_unif)
#     pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
#     pgInterp_unif.interpolate(overwrite=True)
#     x_unif, dataOut_unif = pg.data.select(pgData_unif, z0 = z_xi(z0_coordinate))
#     dataOut_unif = np.squeeze(dataOut_unif)

#     densityFileName_nonunif = str(dataDir+nonunifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_nonunif = pg.GData(densityFileName_nonunif)
#     pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'ms')
#     pgInterp_nonunif.interpolate(overwrite=True)
#     x_nonunif, dataOut_nonunif = pg.data.select(pgData_nonunif, z0 = (z0_coordinate))
#     dataOut_nonunif = np.squeeze(dataOut_nonunif)

#     nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
#     pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
#     pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
#     pgInterp_nonunif_mapc2p.interpolate(overwrite=True)
#     x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pg.data.select(pgData_nonunif_mapc2p, z0 = (z0_coordinate))
#     dataOut_nonunif_mapc2p = np.squeeze(dataOut_nonunif_mapc2p)

#     unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
#     pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
#     pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
#     pgInterp_unif_mapc2p.interpolate(overwrite=True)
#     x_unif_mapc2p, dataOut_unif_mapc2p = pg.data.select(pgData_unif_mapc2p, z0 = z_xi(z0_coordinate))
#     dataOut_unif_mapc2p = np.squeeze(dataOut_unif_mapc2p)

#     unif_Jgeo_filename = str(dataDir+unifFile+'-jacobgeo.gkyl')
#     pgData_unif_Jgeo = pg.GData(unif_Jgeo_filename)
#     pgInterp_unif_Jgeo = pg.GInterpModal(pgData_unif_Jgeo, polyOrder, 'ms')
#     pgInterp_unif_Jgeo.interpolate(overwrite=True)
#     x_unif_Jgeo, dataOut_unif_Jgeo = pg.data.select(pgData_unif_Jgeo, z0 = z_xi(z0_coordinate))
#     dataOut_unif_Jgeo = np.squeeze(dataOut_unif_Jgeo)

#     nonunif_Jgeo_filename = str(dataDir+nonunifFile+'-jacobgeo.gkyl')
#     pgData_nonunif_Jgeo = pg.GData(nonunif_Jgeo_filename)
#     pgInterp_nonunif_Jgeo = pg.GInterpModal(pgData_nonunif_Jgeo, polyOrder, 'ms')
#     pgInterp_nonunif_Jgeo.interpolate(overwrite=True)
#     x_nonunif_Jgeo, dataOut_nonunif_Jgeo = pg.data.select(pgData_nonunif_Jgeo, z0 = (z0_coordinate))
#     dataOut_nonunif_Jgeo = np.squeeze(dataOut_nonunif_Jgeo)

#     dataOut_unif = np.abs(dataOut_unif / dataOut_unif_Jgeo)
#     dataOut_nonunif = np.abs(dataOut_nonunif / dataOut_nonunif_Jgeo)
#     dataOut_diff = np.abs((dataOut_unif - dataOut_nonunif))

#     # Need velocity space grids or to convert to edges of z for plotting
#     fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18,6))

#     norm = LogNorm(vmin = 1e-6, vmax = np.max(dataOut_unif))  # Create a LogNorm instance

#     pcolormesh1 = ax1.pcolormesh(x_unif[1]/vti, x_unif[2]/mui0, dataOut_unif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot
#     # Label the axes
#     ax1.set_ylabel('mu / mui0')
#     ax1.set_xlabel('vpar / vti')
#     ax1.set_title('Frame '+str(frameNum))

#     pcolormesh2 = ax2.pcolormesh(x_nonunif[1]/vti, x_nonunif[2]/mui0, dataOut_nonunif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
#     # Label the axes
#     ax2.set_ylabel('mu / mui0')
#     ax2.set_xlabel('vpar / vti')
    
#     pcolormesh3 = ax3.pcolormesh(x_unif[1]/vti, x_unif[2] / vti, dataOut_diff.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh3, ax=ax3)  # Add a colorbar to the plot
#     # Label the axes
#     ax3.set_ylabel('vpar / vti')
#     ax3.set_xlabel('Absolute difference between distribution functions Z, cylindrical coordinate (m)')

#     ax2.title.set_text('z = '+str(z_xi(z0_coordinate))+' m')

#     locName = np.round(z0_coordinate*100)
#     figName = 'distf_z'+str(locName)+'_'+str(frameNum)
#     if save_figure_as_file:
#       fig.savefig(outDir+figName+figureFileFormat)
#       plt.close()
#       print('Saved figure as '+figName+figureFileFormat)  
#     else:
#       plt.show() 

#     # On axis difference
#     reduce_set_unif = dataOut_unif[:,0]
#     reduce_set_nonunif = dataOut_nonunif[:,0]
#     x_axis = x_unif[1]/vti

#     # I need to interpolate for x 
#     plt.plot(x_axis[1:], reduce_set_unif, label='Uniform grid')
#     plt.plot(x_axis[1:], reduce_set_nonunif, label='Nonuniform grid')
#     plt.plot(x_axis[1:], np.abs(reduce_set_unif - reduce_set_nonunif), 'k--', label='Absolute difference')
#     plt.xlabel('vpar / vti')
#     plt.ylim(np.max(reduce_set_nonunif)* 1e-6, np.max(reduce_set_nonunif)*1.5)
#     plt.xlim(-2, 2)

#     plt.ylabel('f')
#     # log scale y axis
#     plt.yscale('log')
#     plt.legend()
#     plt.title('On axis difference, mu=0, frame '+str(frameNum)+ ', z = '+str(z_xi(z0_coordinate))+' m')

#     figName = 'distf_z'+str(locName)+'mu0'+str(frameNum)
#     if save_figure_as_file:
#       plt.savefig(outDir+figName+figureFileFormat)
#       plt.close()
#       print('Saved figure as '+figName+figureFileFormat)
#     else:
#       plt.show()

#   if plot_distf_at_z:
#     for z0 in z_loctions:
#       plot_distf_at_z_eq(z0)


# # # process_frame(0)
# # # Number of processes to run in parallel
# # num_processes = multiprocessing.cpu_count()
# # print('Number of processes: ', num_processes)

# # # Create a pool of processes
# # pool = multiprocessing.Pool(processes=num_processes)

# # # Map the frame_arr to the pool of processes
# # pool.map(process_frame, frame_arr)

# # # Close the pool to prevent any more tasks from being submitted
# # pool.close()

# # # Wait for all processes to finish
# # pool.join()
  