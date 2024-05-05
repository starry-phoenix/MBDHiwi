"""
This file lists all constants used for the computations. 
The names are (more or less) equivalent to the ones used in the Code by Buffo_et_al. (2018) 
available on GitHub (https://github.com/jbuffo/SlushFund-2.0-interactive).
They use several formulations and constant provided by Griewank and Notz (2013)
"""
import userinput
## Define initial values
S_sw = 34.0 #34.0
T_air = 260.0  # numerical instabilities arise for too low
T_air_Stefan = 265 #300.0 # FREEZING # MELTING
phi_ini = 1 # FREEZING
phi_ini_Stefan = 1  # FREEZING

###
beta = 0.0005836  # Coefficient for density dependence on salinity
kappa = 1.37 * 10 ** (-7)  # Thermal Diffusivity
mu = 1.88 * 10 ** (-3)  # Viscosity
Ra_c = 1.01  # Critical Rayleigh
alpha = 1.56 * 10 ** (-3)  # Linear coeff for Rayleigh number driven advection
Tm_w = 273.15 - 1.853 * S_sw / 28.0    #273.15  # Melt temperature of pure ice [K]
k_i = 2.0  #                            # Thermal conductivity (ice) [W/m/K]
k_br = 0.6  # Thermal conductivity (brine) [W/m/K]
k_w = k_br #2.0
D_s = 2 * 10 ** (-9)  # Diffusivity for Salt
c_br = 3985  # Specific heat of seawater (J/kg/K)
c_i = 2000  # Specific heat of ice
c_w = c_br #4200  # specific heat of water
L = 334774  # Latent heat of fusion ice<->water (J/Kg)
rho_i = 917  # Density of Ice (Kg/m^3)
rho_br = 1028  # Density of Ocean (used in volume averaging - 1D grav. drainage uses delta S) 34ppt NaCl-1027                                     #  12.3ppt MgSO4-1012, 100pppt-1103, 282ppt-1323
rho_w = rho_br #1000
m = 2  # Cementation exponent for Archies equation
g = 9.8  # Earth Gravity
phi_c = 0.06823  # Critical Porosity
P_s = 1 / 100  #    Partition coefficient
a_phi = 0.0000059
b_phi = 1 / a_phi