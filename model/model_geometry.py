import sys
import userinput

def set_model_geometry(geom=2):
    """
    Set up model geometry
    Arguments:
    --------------------
        1      is a test case secnario
        2 is accroding to W3 in Buffo et al. (2018)

    Returns:
    ----------------------
    nz    number of nodes [m]
    dz    cell size       [m]
    Z     total height    [m]
    """
    dz = userinput.dz
    if geom == 1:
        Z = 1
        nc = int(Z / dz)
        nz = int(nc + 1)
    elif geom == 2:
        Z = 1
        nc = int(Z / dz)
        nz = int(nc + 1)
    else:
        print("Requested geometry not available")
    return dz, nz, Z
