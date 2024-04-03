import numpy as np
import math
import userinput
import scipy.optimize as opt
from scipy.special import erfc
import matplotlib.pyplot as plt

if userinput.constants == "real":
    from model.constants_real import T_air_Stefan, Tm_w, L, c_i, k_i, rho_i, c_br, k_br, rho_br, D_s, S_sw
else:
    from model.constants_debug import T_air_Stefan, Tm_w, L, c_i, k_i, rho_i, c_br, k_br, rho_br, D_s, S_sw

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


def analytical_solution_2phase_stefan(depth_stefan, t, dz,nz):
    """
    Analytical solution for the 2-phase Stefan problem.

    Parameters:
    t : time
    x : position
    alpha_l : thermal diffusivity of the liquid phase
    alpha_s : thermal diffusivity of the solid phase
    T_m : melting temperature
    T_i : initial temperature of the solid
    T_0 : boundary temperature at x = 0

    Returns:
    T : temperature at position x and time t
    s : position of the phase change interface at time t
    """
    nz_depth = int(np.absolute(depth_stefan)/dz)
    z = np.linspace(0,np.absolute(depth_stefan), nz_depth)
    T = np.zeros(nz) *Tm_w  #0 
    eta = z/np.sqrt(4*D_s*t)
    T_b = T_air_Stefan  # initial temperature of the solid
    T_inf = Tm_w  # melting temperature of the solid at z=inf in a semi-infinite domain
    eps = lambda k: np.sqrt(D_s/k)
    F_x =lambda x: np.pi**(1/2)*x*np.exp(x**2)*erfc(x)
    C_fi = lambda lamb: S_sw*F_x(lamb)/(1- F_x(lamb))
    S = S_sw + C_fi(depth_stefan)   # C_h
    T_h = Tm_w - 1.853 * S / 28.0
    if eta < depth_stefan:
        # In the solid phase
        T = T_b + (T_h - T_b) * math.erf(eps(k_i)*eta) / math.erf(eps(k_i)*depth_stefan)
    else:
        # In the liquid phase
        T = T_inf + (T_h - T_inf) * erfc(eps(k_br)*eta) / erfc(eps(k_br)*depth_stefan)
        C = S_sw + (S-S_sw) * erfc(eta)/erfc(depth_stefan)

    return T, C

def stefan_problem_twophase(t):
    lamb = 0
    gamma = 0
    C_0 = S_sw
    F_x =lambda x: np.pi**(1/2)*x*np.exp(x**2)*erfc(x)
    G_x = lambda x: np.pi**(1/2)*x*np.exp(x**2)*math.erf(x)
    eps = lambda k: np.sqrt(D_s/k)
    T_b = T_air_Stefan
    T_inf = Tm_w
    T_L = Tm_w - 1.853 * S_sw / 28.0
    T_1 = T_L - T_b
    T_0 = T_inf - T_L
    beta = rho_br*c_br/(rho_i*c_i)
    C_fi = lambda lamb: gamma*C_0*F_x(lamb)/(1- F_x(lamb))
    root_fx_lhs = C_fi(lamb)*(beta/F_x(eps(k_br)*lamb) + 1/G_x(eps(k_i)*lamb))
    root_fx_rhs = T_1/G_x(eps(k_i)*lamb) - beta*T_0/F_x(eps(k_br)*lamb) - L/c_i
    root_fx = root_fx_lhs - root_fx_rhs

    lambda_stefan = opt.newton(root_fx, 1.1)
    stefan_depth = 2*lambda_stefan*np.sqrt(D_s*t)

    return stefan_depth

def plot_stefan_temp_twophase():
    dt = 0.1
    t_pass = 0
    T_arr = []
    C_arr = []
    t_pass_arr = []
    for t in range(userinput.iter_max):
        t_pass = t_pass + dt
        t_pass_arr.append(t_pass)
        depth_stefan = stefan_problem_twophase(t_pass)
        T, C = analytical_solution_2phase_stefan(depth_stefan, t_pass, userinput.dz, userinput.nz)
        T_arr.append(T)
        C_arr.append(C)

    plt.grid()
    plt.plot(t_pass_arr, T_arr, label="Temperature")
    plt.xlabel("Time")
    plt.ylabel("Temperature")
    #plt.plot(t_pass_arr, C_arr, label="Salinity")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    plot_stefan_temp_twophase()