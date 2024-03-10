from model.initial_boundary_conditions import boundary_condition
import numpy as np
import userinput

if userinput.constants == "real":
    from model.constants_real import T_air_Stefan, Tm_w
else:
    from model.constants_debug import T_air_Stefan, Tm_w

def apply_boundary_condition(
    argument,
    X_initial,
    source,
    factor1,
    factor3,
    a,
    Delta_W,
    w,
    nz,
    t_passed,
    S_IC,
    top_temp,
    Stefan,
    Buffo=False, 
    bc_neumann=None,
):
    """
    creates right hand side of matrix equation considering source terms
    Arguments:
    -------
    argument    either 'salinity'' for salt equation or 'temperature' for temperature equation
    X_initial   value of X at last time step
    factor1
    factor3
    Delta_W     difference of ice volume fraction if this and last time step
    nz          number of computational nodes
    t_passed    time passed in seconds

    Returns:
    -------
    B           right hand side of the equation
    """
    B = X_initial
    if argument == "salinity":
        B = B - factor3 * Delta_W
        S_bct, S_bcb = boundary_condition(argument, t_passed, S_IC)
        # as in Buffos code
        B[0] = B[0] + factor1[0] * X_initial[0]
        B[-1] = B[-1] + factor1[-1] * S_bcb  # Dirichlet

    elif argument == "temperature":
        B = B - factor3 * Delta_W  # + 1 / a * source  # check if freezing or melting
        (T_bct, T_bcb) = boundary_condition(argument, t_passed, S_IC, top_temp=top_temp)
        # if Stefan is True:        # for MELTING
        #     B[0] = T_air_Stefan
        #     B[-1] = Tm_w
        if Stefan is True:
            if bc_neumann is not None:  
                B[0] = B[0] - 2*factor1[0]*bc_neumann  
            else:        # for FREEZING
                B[0] = T_air_Stefan #T_air_Stefan #B[-1] + factor1[-1] * T_air_Stefan #T_air_Stefan #B[-1] + factor1[-1] * T_air_Stefan #T_air_Stefan # + factor1[0]*T_air_Stefan
            B[-1] = Tm_w # Dirichlet
        elif Buffo is True:
            if bc_neumann is not None:
                 B[0] = B[0] - 2*factor1[0] * bc_neumann
            else: 
                B[0] = T_air_Stefan
            B[-1] = B[-1] + factor1[-1] * Tm_w
        else:
           # non-pragmatic dirichlet / similar to Buffos code
            B[0] = B[0] + factor1[0] * T_air_Stefan
            B[-1] = B[-1] + factor1[-1] * Tm_w
    else:
        print("argument for BC not available")
    return B

# not relevant for our topic
def correct_for_brine_movement(argument, X_initial, w, t_passed, nz, S_IC, top_temp):
    X_wind = np.zeros(nz)
    sign_w = np.sign(w)
    ones_w = np.ones(nz) * sign_w
    twos_w = np.ones(nz) * 2 * sign_w
    for i in range(1, nz - 1):
        if w[i] != 0:
            w_neg = (
                (ones_w[i] - abs(ones_w[i]))
                / twos_w[i]
                * (X_initial[i + 1] - X_initial[i])
            )  # w < 0
            w_pos = (
                (abs(ones_w[i]) + ones_w[i])
                / twos_w[i]
                * (X_initial[i] - X_initial[i - 1])
            )  # w > 0
            X_wind[i] = w_neg + w_pos
        else:
            X_wind[i] = 0
    if argument == "salinity":
        S_bcb, _ = boundary_condition(argument, t_passed, S_IC)
        if w[0] > 0:
            X_wind[0] = 0  # Neumann BC
        elif w[0] < 0:
            X_wind[0] = X_initial[1] - X_initial[0]
        else:
            X_wind[0] = 0

        if w[-1] > 0:
            X_wind[-1] = X_initial[-1] - X_initial[-2]
        elif w[-1] < 0:
            X_wind[-1] = S_bcb - X_initial[-1]
        else:
            X_wind[-1] = 0

    elif argument == "temperature":
        (T_bct, T_bcb) = boundary_condition(argument, t_passed, S_IC, top_temp=top_temp)
        if w[0] > 0:
            X_wind[0] = X_initial[0] - T_bct
        elif w[0] < 0:
            X_wind[0] = X_initial[1] - X_initial[0]
        else:
            X_wind[0] = 0

        if w[-1] > 0:
            X_wind[-1] = X_initial[-1] - X_initial[-2]
        elif w[-1] < 0:
            X_wind[-1] = T_bcb - X_initial[-1]
        else:
            X_wind[-1] = 0
    else:
        print("argument for BC not available")
    return X_wind
