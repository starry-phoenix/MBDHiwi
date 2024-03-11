"""
Initial conditions, tolerance values of state variables (T,S,P) can be set. 
Geometry type can be selected according to constants_real (geom=2) 
                                        or constants_debug(geom=1) settings.
Maximum iteration number and time step can be varied. (simulation time = iter_max*dt)
"""

constants = "real"          # "real" OR "debug"
if constants=="real":
    from model.constants_real import T_air_Stefan, Tm_w, rho_br, c_br,k_br
    topT_stefan =T_air_Stefan       # Temperature Initial Condition (for T_IC=T_Stefan)
else:
    from model.constants_debug import T_air_Stefan, Tm_w, rho_br, c_br,k_br
    topT_stefan =T_air_Stefan               # Temperature Initial Condition (for T_IC=T_Stefan)

geom = 2            # OR geom = 1 or 2
# Set maximum iteration number, first time step
iter_max = 25000
# set true when comparing model to analytical (Stefan) model
Stefan = True 
Buffo = True
liq_rel = "Normal"    # OR "Normal" OR "FREZCHEM"

dz = 0.01         # grid size
bc_condition = "Dirichlet"  # OR "Neumann" OR "Robin" 
# Define tolerance values for convergence
T_tol = 0.01
S_tol = 0.01
phi_tol = 0.001

# Temperature, Salinity and Partion coeff IC
T_IC='Tm_w'     # OR "T(S)" OR "T271.25" OR 'T_Stefan'
S_IC='S33'           # OR "S_linear" OR "SX" where X is a salinity value; ex S33 with salinity being 33
P_IC='P1'           # OR "P1" OR "P_Stefan" OR 'P0'

# Temperature BC for top layer   
top_temp = "Stefan"         # OR "T_const_250" OR "T_const_260" OR "T_const_265" OR "T_W3" OR "Stefan"      

# for constants = "debug": initial_boundary_conditions.py->L27->replace T_air_Stefan with Tm_w and vice-versa
def fourier_number_timestep():
    alpha = rho_br*c_br
    beta = k_br
    dt = 0.5*alpha*dz**2/beta

    return dt

# Time-step 
dt =47 #0.1*fourier_number_timestep()  # 4.8e-7 for ones and zeros 0.5e-5