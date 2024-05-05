import numpy as np
import re
import userinput

if userinput.constants == "real":
    from model.constants_real import (
    S_sw,
    T_air_Stefan,
    Tm_w
)
else:
    from model.constants_debug import (
    S_sw,
    T_air_Stefan,
    Tm_w
)

def temperature_IC(T_IC, nz):
    Tm = compute_melting_temperature_from_salinity(S_sw)
    Tm_sw = 273.15 - 1.918295  # compute_melting_temperature_from_salinity(S_sw)
    if T_IC == "T(S)":
        T = (
            np.ones(nz, dtype=np.float64) * Tm
        )  #  compute_melting_temperature_from_salinity(S)
    #     T = np.ones(nz, dtype=np.float64) * Tm + 0.001  # slightly above freezing
    elif T_IC == "T_Stefan":
        T = np.ones(nz, dtype=np.float64) * userinput.topT_stefan
    elif T_IC == "T271.25":
        T = np.ones(nz, dtype=np.float64) * 271.25
    elif T_IC=="T250":
        T = np.ones(nz, dtype=np.float64) * 250.0
    elif T_IC=='Tm_w':
        T = np.ones(nz, dtype=np.float64) * Tm_w
    else:
        print("T_IC option not available in initial condition")

    return T

def salinity_IC(S_IC, nz):

    if S_IC =="S_linear":
        S = np.linspace(200, S_sw, nz)
    else:
        S_value = re.findall('[0-9]+$', S_IC)
        if len(S_value)==0:
            raise Exception('S_IC option not available in initial condition')
        else:
            S = np.ones(nz, dtype=np.float64) * float(S_value[0])
    
    return S

def liquidfraction_IC(P_IC, nz):
    #phi = update_liquid_fraction(T, S, phi, H, H_solid, nz, Stefan=False, method="likebuffo")
    if P_IC == "P1":
        phi = np.ones(nz, dtype=np.float64) * 1
    elif P_IC == "P_Stefan":
        phi = np.zeros(nz, dtype=np.float64)
    elif P_IC == 'P0':
        phi = np.zeros(nz, dtype=np.float64)
    else:
        print("P_IC option not available in boundary condition")

    return phi

def set_initial_conditions(Z, nz, T_IC="T0", S_IC="S1", P_IC="P1"):
    """
    Defines initial conditions for: temperature, sea water salinity, porosity#

    Arguments
    ------------------
    nz

    Returns:
    -------------------
    phi       Liquid fraction/porosity profile
    S         Salinity profile
    T         Temperature profile
    w         Brine velocity

    """
    # Brine salinity
    S = salinity_IC(S_IC, nz)
    # Temperature
    T = temperature_IC(T_IC, nz)
    # %% Liquid fraction
    phi = liquidfraction_IC(P_IC, nz)
    # %% brine velocity
    w = np.zeros(nz, dtype=np.float64)

    return T, S, phi, w


def T_W3(t_passed, dt):
    """
    surface temperature changing with time as W3 in Buffo et al.

    Arguments:
    --------
    t_passed    time passe in seconds

    Returns:
    -------
    T_surf      surface temperature at t_passed
    """
    ##Real world boundary condition forcing
    time_temp = np.arange(-2000000, 315576000, dt)
    TtopBC = -9 * np.sin((2 * np.pi / 31557600) * time_temp) - 18.15 + 273.15
    for i in np.range(len(TtopBC)):
        if TtopBC[i] < 250:
            TtopBC[i] = 250
        else:
            TtopBC[i] = TtopBC[i]
    freezedate = (86400 * 149) / dt

    return TtopBC, freezedate

def temperature_BC(t_passed, T_bcb, **kwargs):
        top_temp = kwargs.get("top_temp")
        if top_temp == "T_const_250":
            T_bct = 250.0
        if top_temp == "Stefan":
            T_bct = -1
            T_bcb = 0
        elif top_temp == "T_const_260":
            T_bct = 260.0
        elif top_temp == "T_const_265":
            T_bct = 265.0
        elif top_temp == "T_W3":
            freeze_date = 86400 * 149
            TtopBC = (
                -9 * np.sin((2 * np.pi / 31557600) * ((t_passed + freeze_date)))
                - 18.15
                + 273.15
            )
            if TtopBC < 250:
                TtopBC = 250
            else:
                TtopBC = TtopBC
            T_bct = TtopBC
        else:
            print(
                "Option for top temperature not available. You can define a new option."
            )

        return T_bct, T_bcb

def compute_melting_temperature_from_salinity(S):
    """
    Computes melting temperature at a given salinity
    Arguments
    ---------------
        S   salinity of the brine at melting temperature    [ppt]

    Returns
    ----------
        Tm_b    melting temperature brine/sea water           [K]
    """
    Tm_sw = Tm_w #- 1.853 * S / 28.0  # Tm_w melting temperature pure ice

    return Tm_sw

def boundary_condition(argument, t_passed, S_IC, **kwargs):
    if argument == "temperature":
        
        # temperature at the bottom
        S_value = re.findall('[0-9]+$', S_IC)
        if len(S_value)==0:
            T_bcb = 271.25
        else:
            T_bcb = compute_melting_temperature_from_salinity(float(S_value[0]))

        # temperature at the top
        T_bct, T_bcb = temperature_BC(t_passed, T_bcb, **kwargs)

        return T_bct, T_bcb

    elif argument == "salinity":
        if S_IC =="S_linear":
            S_bcb = 34.0
            S_bct = 100.0
        else:
            S_value = re.findall('[0-9]+$', S_IC)
            if len(S_value)==0:
                raise Exception('S_IC option not available in initial condition')
            else:
                S_bcb = float(S_value[0])
                S_bct = float(S_value[0])

        return S_bct, S_bcb
    
    else:
        print("argument not available. Choose either 'temperature' or 'salinity' ")