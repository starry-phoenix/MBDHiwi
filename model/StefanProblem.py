import numpy as np
import math
import userinput
import scipy.optimize as opt

if userinput.constants == "real":
    from model.constants_real import T_air_Stefan, Tm_w, L, c_i, k_i, rho_i, c_br, k_br, rho_br
else:
    from model.constants_debug import T_air_Stefan, Tm_w, L, c_i, k_i, rho_i, c_br, k_br, rho_br

def stefan_problem(t):
    """
    Computes location of time-dependent phase change interface based on one-phase Stefan problem with Newton iteration

    Argument
    --------------
        t             time for evaluation

    Result
    --------------
        depth_stefan  depth of phase transition interface (anlytically)

    """
    # Stefan number
    #St = c_i * (T_air_Stefan - Tm_w) / L   # for MELTING
    St = c_i * (Tm_w - T_air_Stefan) / L   # for FREEZING
    #St = 1
    # Newton Iteration to obtain lambda
    x = 0
    x_new = 0.1
    while np.absolute(x - x_new) > 0.0003:
        x = x_new
        fx = St * np.exp(-(x**2)) / (np.sqrt(np.pi) * math.erf(x)) - x
        dfxdx = (
            -2 * St * np.exp(-(x**2)) * x / (np.sqrt(np.pi) * math.erf(x))
            - 2 * St * np.exp(-2 * x**2) / (np.pi * math.erf(x) ** 2)
            - 1
        )
        x_new = x - fx / dfxdx
    lam_min = x_new
    # fx = lambda x: St * np.exp(-(x**2)) / (np.sqrt(np.pi) * math.erf(x)) - x
    # dfxdx = lambda x: (
    #         -2 * St * np.exp(-(x**2)) * x / (np.sqrt(np.pi) * math.erf(x))
    #         - 2 * St * np.exp(-2 * x**2) / (np.pi * math.erf(x) ** 2)
    #         - 1
    #     )
    # lam_min = opt.newton(fx, 1.1, fprime=dfxdx)
    alpha = k_i / (c_i * rho_i)  # ice thermal diffusivity
    depth_stefan = 2 * lam_min * np.sqrt(alpha * t)  # positive values
    return depth_stefan


def stefan_temperature(depth_stefan, t, dz, nz):
    """
    Computes temperature distribution of one-phase Stefan problem in the liquid phase
    analytically at all numerical times and all computational nodes

    Arguments
    ----------
        depth_stefan depth of the phase transition (analytical)
        t            evalutation time
        dz           spatially discretization
        nz           number of computational nodes

    Results
    ----------
        T_stefan     Temperature values for all computaional nodes in the liquid phase
    """

    alpha =k_i / (c_i * rho_i)  # ice thermal diffusivity
    nz_depth = int(np.absolute(depth_stefan)/dz)
    z = np.linspace(0,np.absolute(depth_stefan), nz_depth)
    T_stefan = np.ones(nz) *Tm_w  #0 
    if np.absolute(depth_stefan) == nz_depth:
        T_stefan = np.ones(nz) *Tm_w #0 #userinput.topT_stefan 
    else:       # for FREEZING
        for i in range(nz_depth):
            T_stefan[i] = T_air_Stefan - (T_air_Stefan - Tm_w) * (
                math.erf(z[i] / (2 * np.sqrt(alpha * t)))
                / (math.erf(depth_stefan / (2 * np.sqrt(alpha * t))))
            )
    ## for MELTING
    # else:
    #     for i in range(nz_depth):
    #         T_stefan[i] = T_air_Stefan - (T_air_Stefan - Tm_w) * (
    #             math.erf(z[i] / (2 * np.sqrt(alpha * t)))
    #             / (math.erf(depth_stefan / (2 * np.sqrt(alpha * t))))
    #         )
    return T_stefan