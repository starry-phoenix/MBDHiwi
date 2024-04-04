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

S_sw = S_sw/10.0
# TODO: check if the concentration is dependent on previous iterations or vice-versa
# TODO: implement mushy layer consideration with liquid fraction check paper


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


def stefan_temperature_twophase(depth_stefan, t, dz,nz):
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
    T = np.ones(nz)*Tm_w  #0 
    C = np.zeros(nz)  #0
    eta = z/np.sqrt(4*D_s*t)
    T_b = T_air_Stefan  # initial temperature of the solid
    T_inf = Tm_w  # melting temperature of the solid at z=inf in a semi-infinite domain
    eps = lambda k: np.sqrt(D_s/k)
    kappa_i = k_i/(rho_i*c_i)
    kappa_br = k_br/(rho_br*c_br)
    F_x =lambda x: np.pi**(1/2)*x*np.exp(x**2)*erfc(x)
    C_fi = lambda lamb: S_sw*F_x(lamb)/(1- F_x(lamb))
    S = S_sw + C_fi(depth_stefan)   # C_h
    T_h = Tm_w - 1.853 * S / 28.0
    for i in range(nz_depth):
        if eta[i] < depth_stefan:
            # In the solid phase
            T[i] = T_b + (T_h - T_b) * math.erf(eps(kappa_i)*eta[i]) / math.erf(eps(kappa_i)*depth_stefan)
        else:
            # In the liquid phase
            T[i] = T_inf + (T_h - T_inf) * erfc(eps(kappa_br)*eta[i]) / erfc(eps(kappa_br)*depth_stefan)
            C[i] = S_sw + (S-S_sw) * erfc(eta[i])/erfc(depth_stefan)

    return np.array(T), np.array(C)

def stefan_problem_twophase(t):
    gamma = 1.853/28.0
    C_0 = S_sw
    F_x =lambda x: np.pi**(1/2)*x*np.exp(x**2)*erfc(x)
    G_x = lambda x: np.pi**(1/2)*x*np.exp(x**2)*math.erf(x)
    eps = lambda k: np.sqrt(D_s/k)
    T_b = T_air_Stefan
    T_inf = Tm_w
    T_L = Tm_w - gamma * S_sw
    T_1 = T_L - T_b
    T_0 = T_inf - T_L
    kappa_i = k_i/(rho_i*c_i)
    kappa_br = k_br/(rho_br*c_br)
    beta = rho_br*c_br/(rho_i*c_i)
    C_fi = lambda lamb: gamma*C_0*F_x(lamb)/(1- F_x(lamb))
    root_fx_lhs = lambda lamb:C_fi(lamb)*(beta/F_x(eps(kappa_br)*lamb) + 1/G_x(eps(kappa_i)*lamb))
    root_fx_rhs = lambda lamb: T_1/G_x(eps(kappa_i)*lamb) - beta*T_0/F_x(eps(kappa_br)*lamb) - L/c_i
    root_fx = lambda lamb: root_fx_lhs(lamb) - root_fx_rhs(lamb)
    lambda_stefan = opt.newton(root_fx, 0.1, tol=1e-3, maxiter=100)
    stefan_depth = 2*lambda_stefan*np.sqrt(D_s*t)

    return stefan_depth

def plot_stefan_temp_twophase(z_depth=0.5):
    dt= userinput.dt
    t_passed = 0
    T_arr = []
    C_arr = []
    t_pass_arr = []
    depth_stefan_arr = []
    Z = 1
    nc = int(Z / userinput.dz)
    nz = int(nc + 1)
    t = userinput.iter_max*userinput.dt
    for t in range(0, userinput.iter_max):
        t_passed = t_passed + dt
        depth_stefan = stefan_problem_twophase(t_passed)
        T, C = stefan_temperature_twophase(depth_stefan, t, userinput.dz, nz)
        depth_stefan_arr.append(depth_stefan)
        T_arr.append(T)
        C_arr.append(C)
        t_pass_arr.append(t_passed)
    z = int(z_depth*nz)
    T_arr = np.array(T_arr)
    T_z = T_arr[:,z]
    plt.grid()
    plt.plot(np.array(t_pass_arr)/3600, T_z, label="Temperature")
    plt.xlabel("Time in h")
    plt.ylabel("Temperature")
    #plt.plot(t_pass_arr, C_arr, label="Salinity")
    plt.show()
    return T_arr, C_arr

if __name__ == "__main__":
    T, C = plot_stefan_temp_twophase(0.01)