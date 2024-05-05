from numpy import inf
import numpy as np
from scipy.optimize import newton, root, bisect
import userinput

if userinput.constants == "real":
    from model.constants_real import c_i, L, Tm_w
else:
    from model.constants_debug import c_i, L, Tm_w

from decimal import *
getcontext().prec = 4

def compute_melting_temperature_from_salinity(S, liq_rel="Normal"):
    """
    Computes melting temperature at a given salinity
    Arguments
    ---------------
        S   salinity of the brine at melting temperature    [ppt]

    Returns
    ----------
        Tm_b    melting temperature brine/sea water           [K]
    """
    if liq_rel not in ["Normal", "FREZCHEM"]:
        raise TypeError("liquid relation not available")
    if liq_rel == "Normal":
        Tm_sw = Tm_w*np.ones(S.shape) #- 1.853 * S / 28.0  # Tm_w melting temperature pure ice
    elif liq_rel == "FREZCHEM":
        Tm_sw = Tm_w -(-(9.1969758*(1e-05)*S**2)-0.03942059*S+272.63617665)    # FREZCHEM seawater approximation
    else:
        pass

    return Tm_sw

def update_liquid_fraction(
    T, S, phi_k, H, H_solid, nz, Stefan=False, method="likebuffo"
):
    """
    Computes liquid fraction/porosity from enthalpy

    Arguments
    -------------------
      H       enthalpy of the system,
      H_solid     enthalpy of the system if everything was solid ice
      nz      number of computational nodes
    Returns
    -----------------
      phi     liquid fraction [-]
    """

    phi = np.ones(nz, dtype=np.float64)
    alpha = 1.853 / 28.0

    for i in range(nz):
        if Stefan is True:        # for MELTING
            if H[i] <= H_solid[i]:
                phi[i] = 0
            elif H[i] > (H_solid[i] + L):
                phi[i] = 1
            else:
                T_diff = T[i] - compute_melting_temperature_from_salinity(S[i])
                #assert(T[i] == compute_melting_temperature_from_salinity(S[i])), "T - Tm = {}".format(T_diff)  #T[i] = H_solid[i]/c_i
                phi[i] = (H[i] - H_solid[i]) / L
        # if Stefan is True:
        #     if H[i] < (H_solid[i] - L):       # for FREEZING
        #         phi[i] = 0
        #     elif H[i] > (H_solid[i]):
        #         phi[i] = 1
        #     # else:
        #     #     phi[i] = (H[i] - H_solid[i]) / L
        # if Stefan is True:  # melting
        #     H_melt = Tm_w + L
        #     if H[i] > c_i * H_melt:
        #         phi[i] = 1
        #     else:
        #         phi[i] = (H[i] - H_melt) / L
        else:
            if H[i] <=H_solid[i]:
                phi[i] = 0
            elif H[i] > (H_solid[i] + L):
                phi[i] = 1
            else:
                phi[i] = (H[i] - H_solid[i]) / L
            assert (
                phi[i] >= 0 or phi[i] <= 1
            ), "liquid fraction has non physical value"
    return phi * phi_control_for_infinite_values(phi), T


def update_enthalpy_solid_state(S, nz, liq_rel="Normal"):
    """
    compute enthalpy for the cell completely frozen, at a given salinity
    Arguments:
    ----------
      S   salinity [ppt]
      nz  number of computational nodes

    Returns
    ---------
      H_solid   solid state enthalpy
    """
    Tm_b = compute_melting_temperature_from_salinity(S, liq_rel=liq_rel)
    H_solid = np.zeros(nz, dtype=np.float64)
    H_solid = c_i * (Tm_b)
    return H_solid


def update_enthalpy(T, S, phi, nz, method="likebuffo"):
    """
    compute enthalpy of the system

    Arguments
    ---------
      T       temperature [K]
      S       salinity [ppt]
      phi     porosity [-]
      nz      number of computational nodes

    Returns
    --------
      H       enthalpy of the system
    """
    if method not in ["likebuffo"]:
        raise TypeError("method for enthalpy not available")
    if method == "likebuffo":
        H =  L * phi + c_i * T 
    else:
        pass

    return H


def phi_control_for_infinite_values(phi):
    """
    for the some computations, division by phi is required
    to circumvent infinite values for low liquid fractions

    Arguments
    -----------
      phi liquid fraction

    Returns
    -----------
      phi_control vector to multiply with phi
    """
    one_test = np.ones_like(phi, dtype=np.float64)
    phi_not_low = np.zeros_like(phi, dtype=np.float64)
    for i in range(len(phi)):
        if phi[i] > 0.000001:
            phi_not_low[i] = one_test[i] / phi[i]
        else:
            phi_not_low[i] = inf
    phi_control = np.ones_like(phi, dtype=np.float64)
    phi_control[np.isneginf(phi_not_low)] = 0
    phi_control[np.isposinf(phi_not_low)] = 0
    return phi_control


#### NEWTON ITERATION
def H_function(x, Tk, H_Sl):
    a = 0.0000059
    b = 1 / a
    return (
        Tk * c_i + L * 0.5 * (np.tanh(a * x - (H_Sl + L / 2) / b) + 1) - x
    )  #  np.piecewise(x, [x<H_Sl, H_Sl+L>=x>= H_Sl, x>H_Sl+L], [0, lambda x : (x-H_Sl)/L, 1]))-x


def H_function_derivate(x, Tk, H_Sl):
    a = 0.0000059
    b = 1 / a
    return (0.5 * a * L) / (np.cosh(a * x - (L / 2 + H_Sl) / b) ** 2) - 1
    # (np.piecewise(x, [x<H_Sl, H_Sl+L>=x>= H_Sl, x>H_Sl+L], [-1, 0 , -1]))


def H_newton_iteration(Tk, H_Sl):
    """
    newton iteration for each pair of (Tk,Sl) to find enthalpy
    Arguments
    -------------------
    Tk:    temperature scalar
    Sl:    brine salinity scalar

    Returns
    ------------
    Hkl:  enthalpy of the system at position k,l
    """
    x = 1
    Hkl = newton(
        H_function, x0=x, fprime=H_function_derivate, args=(Tk, H_Sl), maxiter=100
    )
    return Hkl


def phi_func(Hkl, H_Sl):
    """
    function to compute phi
    Arguments:
    ---------------
    Hkl:  enthalpy at k,l
    Sl:   salinity at l
    """
    # phi = 0.5*(np.tanh(a_phi*Hkl-(H_Sl+L/2)/b_phi)+1)
    # phi = np.piecewise(Hkl, [Hkl<H_Sl, H_Sl+L>=Hkl>= H_Sl, Hkl>H_Sl+L],\
    # [0, lambda Hkl : (Hkl-H_Sl)/L, 1])
    a = 0.0000059
    b = 1 / a
    phi = 0.5 * (np.tanh(a * Hkl - (H_Sl + L / 2) / b) + 1)
    return phi


def H_update(T, H_solid, nz):
    """
    computes enthalpy from T and S
    """
    H = np.zeros(nz)
    for k in range(nz):
        Tk = float(T[k])
        # phik = float ( phi[k] )

        H_Sk = H_solid[k]
        Hk = H_newton_iteration(Tk, H_Sk)

        H[k] = Hk
    return H


def H_newton_manual(T_k, H_Sl):
    # Newton Iteration to obtain lambda
    x = T_k * c_i + L
    x_new = x + 1000
    while np.absolute(x - x_new) > 0.0003:
        x = x_new
        fx = (
            T_k * c_i
            + L
            * np.piecewise(
                x,
                [x < H_Sl, H_Sl + L >= x >= H_Sl, x > H_Sl + L],
                [0, lambda x: (x - H_Sl) / L, 1],
            )
            - x
        )
        dfxdx = np.piecewise(
            x, [x < H_Sl, H_Sl + L >= x >= H_Sl, x > H_Sl + L], [-1, 0, -1]
        )
        x_new = x - fx / dfxdx
    return x_new


def H_bisect_search(Tk, H_Sl):
    x1 = Tk * c_i + L * 0
    x2 = Tk * c_i + L * 1
    try:
        root = bisect(H_function, a=x1, b=x2, args=(Tk, H_Sl), xtol=0.1)
    except RuntimeError:
        low = H_function(x1, Tk, H_Sl)
        high = H_function(x2, Tk, H_Sl)
        if low < high:  # low is closer to zero
            root = x1
        elif high < low:  # high is closer to zero
            root = x2
    print(root)
    return root
