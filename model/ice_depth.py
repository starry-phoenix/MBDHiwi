import numpy as np

def locate_ice_ocean_interface(phi, dz, nz, **kwargs):
    """
    Locate ice ocean interface, based on liquid fraction
    equivalent ice thickness

    Arguments
    -------------
        phi             liquid fraction [-]
        dz              spatial discreitization [m]
        nz              number of computational nodes
        ***kwargs       Validation with Stefan problem: Stefan = True

    Results
    -----------
        if_depth        location of the ice-water interface/sea ice total thickness [m]
        if_depth_index  index of the 'transition cell' from ice to ocean (freezing) or water to ice (Melting)
    """
    # initialize variable for for-loop

    Stefan = kwargs.get("Stefan", True)
    Z = (nz - 1) * dz
    depth = np.linspace(dz, Z + dz, nz)
    ice_depth_index = 0
    ice_depth = 0
    if Stefan is False:  # Freezing
        for i in range(nz):
            if phi[i] < 1:  # if cell is not completely liquid it counts as sea ice
                ice_depth = depth[i]
                ice_depth_index = i + 1
            else:
                ice_depth = ice_depth
                ice_depth_index = ice_depth_index
    elif Stefan is True:  # Melting water thickness
        for i in range(nz):
            if phi[i] >= 1-0.05:   # if phi[i] >= 1-0.01 for FREEZING and if phi[i] <= 0.01
                break
            else:
                ice_depth = depth[i]
                ice_depth_index = i + 1  # number of nodes immersed by water
    return ice_depth, ice_depth_index
