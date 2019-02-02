import numpy as np
import h5py
import matplotlib.pyplot as pl
import matplotlib.colors as mc
from scipy.interpolate import NearestNDInterpolator
from numpy import linspace, logspace

# Define functions, which we need later on. 
# Estimate x, y and z-axis.
def axis_def(content):
    content = content.transpose()
    x_axis = content[0] 
    y_axis = content[1] 
    z_axis = content[2] 
    content = content.transpose()
    return x_axis, y_axis, z_axis
    
# Interpolate x, y, z. 
# nl == number of interpolations in x and y direction
# nh == number of interpolations in z direction
def grid(nl, nh):
    xgrid = linspace(-width, width, nl) 
    ygrid = linspace(-width, width, nl)
    zgrid = linspace(-width_z, width_z, nh)
    XGRID, YGRID, ZGRID = np.meshgrid(xgrid, ygrid, zgrid)
    return XGRID, YGRID, ZGRID, xgrid, ygrid
    
# Important variables and constants, which we need later on.
halonumber   = 208 # Halo 208 of over 9000 in my Simulation
color_range  = 256 # Range of 256 colours for the Plot
contour_lev  = [1, 10, 100]
width    = 0.5 # X and y side width
width_z  = 0.01 # Small z axis depth, we just want a slice at the center
cbar_tic = logspace(-2, 3, 6) #Tics of the colourbar
cbar_pos = [0.8, 0.13, 0.03, 0.77] #Position of the colourbar
cbar_min = -2.0 #Minimum of the colourbar
cbar_max = 3.5 #Maximum of the colourbar
mu = 1.22 #constant to convert in physical units
mp = 1.672623e-24 #constant to convert in physical units

# Take important values and arrays from the file.
ffile ='haloes_013.'+str(halonumber)+'.hdf5'
f     = h5py.File(ffile, 'r')
h           = f['Header'].attrs['HubbleParam']
z           = f['Header'].attrs['Redshift']
l_unit      = f['Header'].attrs['UnitLength_in_cm']
m_unit      = f['Header'].attrs['UnitMass_in_g']
v_unit      = f['Header'].attrs['UnitVelocity_in_cm_per_s']
velo        = np.array(f['PartType0']['Velocities']) #Velocity
dens        = np.array(f['PartType0']['Density']) #Density
cord        = np.array(f['PartType0']['Coordinates']) #Coordinates
halo_cord   = np.array(f['Header']['GroupPos'])
halo_velo   = np.array(f['Header']['GroupVel'])

# Convert the data from code units to physical units.
velo_prev = v_unit*1.0/(np.sqrt(z+1.0))
velo=np.multiply(velo_prev,velo)
velo_prev2 = v_unit*(z+1.0)
halo_velo = np.multiply(velo_prev2,halo_velo)
dens_prev = (m_unit / l_unit**3) * ((z+1.0)**3 * h**2)
dens = np.multiply(dens_prev,dens)
dens_to_ndens = 1.0/(mp * mu)
ndens = np.multiply(dens_to_ndens,dens) 
mass_prev = 1.0e10/h

# Here the 'real' work starts.
# Substract halo coordinate and halo velocity to center the halo.
cord = cord - halo_cord
velo = velo - halo_velo
#Estimate the x, y and z axis for the coordinates and velocities.
cx, cy, cz = axis_def(cord)
vx, vy, vz = axis_def(velo)

# Density and velocity interpolation.
points = np.vstack((cx, cy, cz)).T   
ndensf = NearestNDInterpolator(points, ndens)
vxf = NearestNDInterpolator(points, vx)
vyf = NearestNDInterpolator(points, vy)
vzf = NearestNDInterpolator(points, vz)
# Interpolations for the coordinates(101, 51) -> 101x101x51 cells
# Interpolations for the velocityarrows(22, 1) -> 22x22 arrows
XGRID, YGRID, ZGRID, xgrid, ygrid = grid(101, 51) 
XGRID2, YGRID2, ZGRID2, xgrid2, ygrid2 = grid(22, 1)
# Interpolate velocity in 2D-Array.
VX = vxf(XGRID2,YGRID2,ZGRID2).sum(axis=2) 
VY = vyf(XGRID2,YGRID2,ZGRID2).sum(axis=2)
VZ = vzf(XGRID2,YGRID2,ZGRID2).sum(axis=2)
# Interpolate density in 2D-Array.
NDENS = ndensf(XGRID,YGRID,ZGRID)
NDENS_slice = NDENS.sum(axis=2)/float(51)

# Here the plot commands start.
# First I calculate the size and plot the previously estimated densities
# of the respective cells.
ax=pl.subplot(111)
extent = (xgrid.min(),xgrid.max(),ygrid.min(),ygrid.max())
color_levels = logspace(cbar_min,cbar_max,color_range)
norm = mc.BoundaryNorm(color_levels, color_range)
im = ax.imshow(NDENS_slice, #
               origin='lower',
               norm=norm,
               extent=extent, 
               cmap='jet'
               )
               
# Insert the conturlines. 
# These symobolize the transition to denser regions.
im2 = ax.contour(xgrid,ygrid,NDENS_slice,
                 levels=contour_lev,
                 origin='lower',
                 colors=['black', 'red', 'blue'],
                 linewidths=(1.6, 1.6)
                 )
pl.clabel(im2, inline=6, fontsize=11,fmt='%1.0i')

# Insert arrows into plot.
# These symbolize the velocities of the respective gascell.  
quivers = ax.quiver(XGRID2, YGRID2, VX, VY, color='grey')                
qk = ax.quiverkey(quivers, 0.4, 0.15, 500000, 
               r'$5\, \rm{km}\,/\rm{s}$', 
               labelpos='W',
               coordinates='figure',
               color='black',
               fontproperties={'weight': 'bold', 'size' : '12'}
               )
qk.text.set_backgroundcolor('w')

# Insert the colorbar and labels.
cbar_ax = pl.gcf().add_axes(cbar_pos)
cbar = pl.colorbar(im, ticks=cbar_tic, format='%.0e', cax=cbar_ax)   
cbar.set_label('number density [cm$^{-3}$]', fontsize=15)
ax.set_xlabel(r'$x\, [\rm{ckpc/h}]$', fontsize=15)
ax.set_ylabel(r'$y\, [\rm{ckpc/h}]$', fontsize=15)

# Save the image in a PDF format.
'''
#pl.savefig('Halo_Nr'+str(halonumber)+'.pdf', bbox_inches='tight')
'''
pl.show()