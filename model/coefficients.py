import userinput
from model.liquid_fraction import phi_control_for_infinite_values
import numpy as np

if userinput.constants == "real":
    from model.constants_real import (rho_br,
    rho_i, c_br, c_i, k_i, k_br, P_s,
    L, D_s, m, rho_w, c_w, k_w)
else:
    from model.constants_debug import (rho_br,
    rho_i, c_br, c_i, k_i, k_br, P_s,
    L, D_s, m, rho_w, c_w, k_w)

def update_coefficients(argument, X_initial, w, phi, nz, S_IC):
    """
    Updates of coefficients required to solve the Advection Reaction Diffusion Equation 
    for each time step for temperature or salinity

    Arguments
    --------------
        argument    either 'temperature' or 'salinity'
        X_initial   initial value for salinity or temperature   [ppt] or [K]
        phi         liquid fraction [-]
        w           brine velocity [ms-1]
        nz          number of computational nodes
    
    Returns
    ------------
        a   'temperature': heat capacity | 'salinity': liquid fraction
        c   'temperature': thermal conductivity | 'salinity': salt diffusivity
        d   'temperature': latent heat | 'salinity': factor to determine salinity increase due to liquid fraction decrease

    """
    a = np.zeros(nz)
    b = np.zeros(nz)
    c = np.zeros(nz)
    d = np.ones(nz)
    assert argument in [
        "temperature",
        "salinity",
    ], "invalid input for argument in `update coefficients`"
    if argument == "temperature":
        T = X_initial
        d = d * (rho_i * L)
        if S_IC == "S0" or S_IC == "S1":  # no brine
            a = np.ones(nz, dtype=np.float64) * effective_parameter(
                rho_w * c_w, rho_i * c_i, phi, nz
            )
            b = np.ones(nz, dtype=np.float64) * rho_w * c_w * w
            c = np.ones(nz, dtype=np.float64) * effective_parameter(k_w, k_i, phi, nz)
        else:
            a = np.ones(nz, dtype=np.float64) * effective_parameter(
                rho_br * c_br, rho_i * c_i, phi, nz
            )
            b = np.ones(nz, dtype=np.float64) * rho_br * c_br * w
            c = np.ones(nz, dtype=np.float64) * effective_parameter(k_br, k_i, phi, nz)
    elif argument == "salinity":
        S = X_initial
        a = np.ones(nz, dtype=np.float64) * phi
        b = np.ones(nz, dtype=np.float64) * w
        c = np.ones(nz, dtype=np.float64) * salt_diffusivity(phi, nz)
        d = np.ones(nz, dtype=np.float64) * S * rho_i / rho_br #* (1 - P_s)
    else:
        raise ValueError("Option in coeffients not available")
    return a, b, c, d


def salt_diffusivity(phi, nz):
    """
    Computes effective diffusion coefficient of salt based on Archies law

    Arguments
    --------------------------
        phi     ice volume fraction [-]
        nz      number of computational nodes

    Returns
    ---------------
        D_eff   D_eff effective diffusion coefficient [s2m-1]
    """
    D_eff = np.zeros(nz, dtype=np.float64)
    for i in range(nz):
        D_eff[i] = D_s * phi[i] ** m  # Archies law
    return D_eff


def effective_parameter(A_br, A_i, phi, nz):
    """
    Equation for effective parameters based on brine and ice physical properties based on liquid fraction 
    Arguments
    ---------------------
        A_br    brine property
        A_ i    ice property
        phi     liquid fraction
        nz      number of computational nodes

    Returns:
    ----------------------
        eff     effective coefficient
    """
    eff = np.zeros(nz, dtype=np.float64)
    eff = phi * A_br + (1 - phi) * A_i
    return eff

