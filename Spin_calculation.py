import numpy as np
import h5py

#Some functions we need later on:
#Function to calculate the x, y and z coordinate and the distance
#of each particle in respect to the center.                 
def coordinates(cord): 
    cord_rel_center = cord - halo_pos #substract the halo velocity
    cord_rel_center = cord_rel_center.transpose() 
    cord_x = cord_rel_center[0]  #x-coordinate
    cord_y = cord_rel_center[1]  #y-coordinate
    cord_z = cord_rel_center[2]  #z-coordinate
    cord_rel_center = cord_rel_center.transpose()
    part_dist = np.sqrt((cord_x**2)+(cord_y**2)+(cord_z**2)) 
    return part_dist, cord_x, cord_y, cord_z

#Function to calculate the x, y and z velocity of each particle.                 
def velocities(velo): 
    velo_rel_center = velo - halo_vel 
    velo_rel_center = velo_rel_center.transpose()
    velo_x = velo_rel_center[0]
    velo_y = velo_rel_center[1]
    velo_z = velo_rel_center[2]
    velo_rel_center = velo_rel_center.transpose()
    return velo_x, velo_y, velo_z

#Create an one array that contains the gas and dark matter conents 
#(coordinates or velocity). By doing that I 'merge' both together.
def merging(content, contentDM):
    merged_list  = []
    merged_list.extend(content)
    merged_list.extend(contentDM)
    contentTOT     = merged_list
    contentTOT     = np.array(contentTOT)
    return contentTOT

# Function to calculate the angular momentum which is described by the 
# Length of the angular momentum in x, y and z direction.
def angular_momentum(mass, cord_x, cord_y, cord_z, velo_x, velo_y, velo_z):
    ang_mom_x = np.sum(np.multiply(mass, (cord_y*velo_z) - (cord_z*velo_y)))
    ang_mom_y = np.sum(np.multiply(mass, (cord_z*velo_x) - (cord_x*velo_z)))
    ang_mom_z = np.sum(np.multiply(mass, (cord_x*velo_y) - (cord_y*velo_x)))
    ang_mom   = np.sqrt(ang_mom_x**2 + ang_mom_y**2 + ang_mom_z**2)
    return ang_mom

# Function to calculate the spin. This is a unitless dimension between 
# 0 and 1 that indicates how fast a halo rotates.
def Spin(angul_mom,halo_masses,circ_velo,max_dist):
    Spin = angul_mom / (np.sqrt(2) * halo_masses * circ_velo * max_dist)
    return Spin

# Some constants we need later on: 
mu  = 1.22
mp  = 1.672623e-24
kb  = 1.380658e-16
G   = 6.67408e-8 # Gravitational constant in cm^3/(g*s^2)
M   = 1.9989e33  # Sunmass in gramm
gamma = 5.0/3.0
mass_ofDM   = 6.72669027e-9 # Resolution of dark matter content
halonumber  = 9021 # Our Simulation contains over 9000 haloes
part_limit  = 50 # Halo must contain at least 50 particles to be relevant

# Some arrays we need later on to save our data:
spin_valueTOT 	= np.zeros(halonumber) # Each index represents a halo
spin_valueGAS	= np.zeros(halonumber)
spin_valueDM 	= np.zeros(halonumber)

# Here I loop over all haloes in the Simulation.
for i in xrange(0,halonumber):
    # Take important values and arrays from the file
    try:
	    file='haloes_013.'+str(i)+'.hdf5' 
	    f = h5py.File(file, 'r')
    except:
	    #print 'file does not exist!'
	    continue

    h = f['Header'].attrs['HubbleParam']
    z = f['Header'].attrs['Redshift']
    l_unit  = f['Header'].attrs['UnitLength_in_cm']
    m_unit  = f['Header'].attrs['UnitMass_in_g']
    v_unit  = f['Header'].attrs['UnitVelocity_in_cm_per_s']
    cordGAS = np.array(f['PartType0']['Coordinates'])   # Coordinates
    veloGAS = np.array(f['PartType0']['Velocities'])    # Velocity
    massGAS = np.array(f['PartType0']['Masses'])        # Density
    cordDM  = np.array(f['PartType1']['Coordinates'])   
    veloDM  = np.array(f['PartType1']['Velocities'])     
    massDM  = np.ones(len(cordDM))  # Dark matter particles have all
    massDM  = massDM * mass_ofDM    # the same mass.
    halo_pos= np.array(f['Header']['GroupPos'])     # Position of the halo
    halo_vel= np.array(f['Header']['GroupVel'])     # Velocity of the halo
    virrad  = np.array(f['Header']['Group_R_Crit200'])  # Virialradius
    f.close()

    # We need to convert the data from code units to physical units
    mass_prev   = 1.0e10 / h
    massGAS     = np.multiply(mass_prev, massGAS)
    massDM      = np.multiply(mass_prev, massDM)
    velo_prev   = v_unit / (np.sqrt(z+1.0))
    veloGAS     = np.multiply(velo_prev, veloGAS)
    veloDM      = np.multiply(velo_prev, veloDM)
    velo_prev2  = v_unit * (z+1.0) 
    halo_vel    = np.multiply(velo_prev2, halo_vel)
    
    # Here the 'real' work starts.
    # Determine coordinates of the gas & dark matter halo and calc. 
    # their distance to the center
    part_distGAS, cord_xGAS, cord_yGAS, cord_zGAS = coordinates(cordGAS)
    part_distDM, cord_xDM, cord_yDM, cord_zDM = coordinates(cordDM)
    
    # Only choose particles inside the virial radius
    ids_radius_selectGAS = part_distGAS <= virrad   # True/False array 
    ids_radius_selectDM = part_distDM <= virrad     # True/False array 
    cordGAS = cordGAS[ids_radius_selectGAS]         # Select only 'True'
    massGAS = massGAS[ids_radius_selectGAS]
    veloGAS = veloGAS[ids_radius_selectGAS]
    cordDM  = cordDM[ids_radius_selectDM]
    veloDM  = veloDM[ids_radius_selectDM]
    massDM  = massDM[ids_radius_selectDM]
    
    # To receive an statistical value we need at least part_limit = 50.
    # particles for gas and dark matter. If thats not the case, we calc.
    # the next halo.
    if (len(cordGAS) & len(cordDM)) < part_limit: 
	continue  

    # Determine coordinates inside the virial radius for gas & dark matter.
    part_distGAS, cord_xGAS, cord_yGAS, cord_zGAS = coordinates(cordGAS)
    part_distDM, cord_xDM, cord_yDM, cord_zDM = coordinates(cordDM)
    max_dist = np.max(part_distGAS)
    max_dist_cgs = max_dist*l_unit/((z+1.0)*h) 
    velo_xGAS, velo_yGAS, velo_zGAS = velocities(veloGAS)
    velo_xDM, velo_yDM, velo_zDM = velocities(veloDM)
    
    # Calculate the angular momentum of the gas & dark matter halo.
    ang_momGAS = angular_momentum(massGAS, 
                                  cord_xGAS,
                                  cord_yGAS, 
                                  cord_zGAS, 
                                  velo_xGAS, 
                                  velo_yGAS, 
                                  velo_zGAS
                                  )
    ang_momDM  = angular_momentum(massDM, 
                                  cord_xDM, 
                                  cord_yDM, 
                                  cord_zDM, 
                                  velo_xDM, 
                                  velo_yDM, 
                                  velo_zDM
                                  )
    
    # Merge gas und dark matter halo to a 'total' halo together
    massTOT = merging(massGAS, massDM)
    cordTOT = merging(cordGAS, cordDM)
    veloTOT = merging(veloGAS, veloDM)
    
    # Determine coordinates inside the virial radius for the merged halo.
    # and calculate the angular momentum
    part_distTOT, cord_xTOT, cord_yTOT, cord_zTOT = coordinates(cordTOT)
    velo_xTOT, velo_yTOT, velo_zTOT = velocities(veloTOT)
    halo_massesTOT  = np.sum(massTOT) # Entire mass ist sum of all particles
    halo_massesDM   = np.sum(massDM) 
    halo_massesGAS  = np.sum(massGAS)
    ang_momTOT = angular_momentum(massTOT, 
                                  cord_xTOT, 
                                  cord_yTOT, 
                                  cord_zTOT, 
                                  velo_xTOT, 
                                  velo_yTOT, 
                                  velo_zTOT
                                  )
                                  

    
    # Calculate the circular velocity (always of the total halo) and with 
    # that the spin of gas, dark matter, and total matter.
    halo_mass_cgsTOT = halo_massesTOT * M #in cgs units
    circ_velo        = np.sqrt((G*halo_mass_cgsTOT) / max_dist_cgs)    
    spinTOT = Spin(ang_momTOT, halo_massesTOT, circ_velo, max_dist)
    spinGAS = Spin(ang_momGAS, halo_massesGAS, circ_velo, max_dist)
    spinDM  = Spin(ang_momDM, halo_massesDM, circ_velo, max_dist)
    spin_valueTOT[i] = spinTOT
    spin_valueGAS[i] = spinGAS
    spin_valueDM[i] = spinDM
    
    # Print comands of the spin.
    print 'Spin of total halo nr.'+str(i)+': ', spinTOT
    print 'Spin of gashalo nr.'+str(i)+': ', spinGAS
    print 'Spin of dark matter halo nr.'+str(i)+': ', spinDM    
    
# At the end we can save our data.  
'''
np.savetxt('SpinTOT.out', spin_valueTOT, delimiter=',')                  
np.savetxt('SpinGAS.out', spin_valueGAS, delimiter=',')                  
np.savetxt('SpinDM', spin_valueDM, delimiter=',')                  
'''