"""
This file lists all constants used for the computations. 
The names are (more or less) equivalent to the ones used in the Code by Buffo_et_al. (2018) 
available on GitHub (https://github.com/jbuffo/SlushFund-2.0-interactive).
They use several formulations and constant provided by Griewank and Notz (2013)
"""
import numpy as np
## Define initial values 
S_sw = 0
T_air = 1 #numerical instabilities arise for too low 
T_air_Stefan = -1
phi_ini = 0    # phi_ini = 0 for FREEZING  
phi_ini_Stefan = 0.0

###
beta = 1 #0.0005836                      # Coefficient for density dependence on salinity
kappa = 1 #1.37*10**(-7)                 # Thermal Diffusivity
mu = 1 #1.88*10**(-3)                    # Viscosity
Ra_c = 0 #1.01                           # Critical Rayleigh 
alpha = 1 #1.56*10**(-3)                 # Linear coeff for Rayleigh number driven advection
Tm_w = 0                       # Melt temperature of pure ice [K]
k_i = 1 #2.0  #                            # Thermal conductivity (ice) [W/m/K]
k_br = 1 # 0.6                            # Thermal conductivity (brine) [W/m/K]
k_w =1 # 0.6
D_s = 0 #2*10**(-9)                      # Diffusivity for Salt 
c_br = 1 #3985                           # Specific heat of seawater (J/kg/K)
c_i = 1# 2000                            # Specific heat of ice
c_w =1 #4200                             # specific heat of water
L = 1 #334774                            # Latent heat of fusion ice<->water (J/Kg)
rho_i = 1 #917                           # Density of Ice (Kg/m^3)
rho_br = 1 #1028                        # Density of Ocean (used in volume averaging - 1D grav. drainage uses delta S) 34ppt NaCl-1027                                     #  12.3ppt MgSO4-1012, 100pppt-1103, 282ppt-1323
rho_w = 1 #1000
m = 1 #2                                 # Cementation exponent for Archies equation
g = 0 #9.8                               # Earth Gravity
phi_c =1 # 0.05                          # Critical Porosity
P_s = 1  #-1/100                             # Partition coefficient
a_phi = 1 #0.0000059
b_phi = 1 #1/a_phi

# def ice_density(T):
#   rho_i = 917-0.1403*T
#   return rho_i
# def brine_density(S):
#     rho_br = 1000+0.8*S
#     return rho_br

# def brine_thermalconductivity(T):
#     # Yen et al. 1991, Sea Ice Book
#     k_br = 0.4184 *(1.25+0.03*T +0.00014*T**2)
#     return k_br

# def ice_thermalconductivity(T):
#     # Yen et al. 1991, Sea Ice Book
#     k_i =1.16 *(1.91-8.66*10**(-3)*K+2.97*10**(-5)*T**(-2))
#     return k_i