import userinput
import os
import sys
import time
import matplotlib.gridspec as gridspec
import numpy as np
import math
import pandas as pd
from IPython.display import display, clear_output
import matplotlib.pyplot as plt
import seaborn as sns
from StefanProblem import stefan_problem, stefan_temperature, stefan_problem_twophase, stefan_temperature_twophase
from model.ice_depth import locate_ice_ocean_interface
from model.store import create_directory, set_up_matrices, store_results
from model.model_geometry import set_model_geometry
from model.set_time import set_up_iter, t_total
from model.initial_boundary_conditions import set_initial_conditions, boundary_condition
from model.liquid_fraction import update_enthalpy, update_enthalpy_solid_state, update_liquid_fraction
from model.statevariables_errors import reset_error_for_while_loop, compute_error_for_convergence, reset_error_for_while_loop, initialize_statevariables, overwrite_statevariables, set_statevariables, define_previous_statevariable
from model.solve_T_S import AdvectionDiffusion
import time
from calculate_error import error_norms

sys.path.append('Users\sneha\RWTHMSc\Sem2\MBD_Parttime\GitClone_empty\sea-ice-model-buffo\model\info.txt')

#  Set up the model geometry
if userinput.constants == "real":
    from model.constants_real import T_air_Stefan, Tm_w, L, c_i, k_i, rho_i, phi_c
else:
    from model.constants_debug import T_air_Stefan, Tm_w, L, c_i, k_i, rho_i, phi_c

# [x]NOTE: Added Neumann for Buffo Code and compares well with Matlab code

class SeaIceModel(error_norms):
    def __init__(self):
        geom = userinput.geom
        [self.dz, self.nz, self.Z] = set_model_geometry(geom)
        # Set maximum iteration number, first time step
        self.dt = userinput.dt
        self.iter_max = userinput.iter_max
        self.bc_condition = userinput.bc_condition

        # Define tolerance values for convergence iteration
        self.T_tol = userinput.T_tol
        self.S_tol = userinput.S_tol
        self.phi_tol = userinput.phi_tol

        #dt= 0.5e-5 # 4.8e-7 should be less than tnp1 condition. for the case of unit parameter testing, it should be less than 0.5e-4
        T_IC= userinput.T_IC
        self.S_IC= userinput.S_IC
        P_IC= userinput.P_IC
        self.Stefan = userinput.Stefan
        self.Buffo = userinput.Buffo
        self.liq_rel = userinput.liq_rel
        self.phi_c = phi_c 
        self.critical_depth = 0.01
        if self.bc_condition == "Neumann":
            self.temp_grad = self.dz*(Tm_w - T_air_Stefan)/self.critical_depth     # T_bottom - T_top / depth
        else:
            self.temp_grad = None

        # Initialize
        self.error_temperature, self.error_temperature_sum, self.error_temperature_sum_weighted = np.zeros(self.nz), np.zeros(self.iter_max), np.zeros(self.iter_max)
        self.thickness_index_total, depth_stefan_all = np.zeros(self.iter_max), np.zeros(self.iter_max)
        self.T_Stefan_diff, self.T_Stefan_prev, self.T_k_diff, self.T_k_prev = [], np.zeros(self.nz), [], np.zeros(self.nz)
        self.T_Stefan_list, self.T_k_list, self.T_k_buffo_list,self.thickness_list, self.thickness_list_Buffo = [], [], [], [], []
        
        [self.iter_max, self.dt, self.t_passed] = set_up_iter(self.iter_max)

        [self.T, self.S, self.phi, self.w] = set_initial_conditions(self.Z, self.nz, T_IC, self.S_IC, P_IC)  # define initial conditons for temperature, brine salinity and liquid fraction
        self.H_solid = update_enthalpy_solid_state(self.S, self.nz, liq_rel=self.liq_rel)  # compute solid enthalpy
        self.H = update_enthalpy(self.T, self.S, self.phi, self.nz)
        thickness, thickness_index = locate_ice_ocean_interface(self.phi, self.dz, self.nz, Stefan= self.Stefan)
        [self.all_T, self.all_S, self.all_phi, self.all_H, self.all_H_solid, self.all_w, self.all_thick, self.all_t_passed]= \
            set_up_matrices(self.iter_max, self.nz)
        [self.all_T, self.all_S, self.all_phi, 
         self.all_H, self.all_H_solid, self.all_w, 
         self.all_thick, self.all_t_passed] = \
                store_results(self.T, self.S, self.phi, self.H, self.H_solid, 
                              self.w, thickness, self.t_passed, self.all_T, 
                              self.all_S, self.all_phi, self.all_H, self.all_H_solid, 
                              self.all_w, self.all_thick,self.all_t_passed, 0
                        ) 
        self.folder_name = "Temperature_" + self.bc_condition + "_"+str(self.dz) + "_" + str(self.dt) + "_" + str(self.iter_max) 
        if not os.path.isdir(self.folder_name):
            os.makedirs(self.folder_name)

    def T_running(self,fig, ax1, T_Stefan, T_k, T_k_buffo=None, count=0): # script to visualize analytical and numerical solution updates
        ax1.plot(T_k,np.linspace(0,-self.Z, self.nz), 'r--')
        ax1.plot( T_Stefan,np.linspace(0,-self.Z, self.nz), 'k')
        if T_k_buffo is not None:
            ax1.plot(T_k_buffo,np.linspace(0,-self.Z, self.nz), 'b-.',alpha=0.5)
        else:
            pass
        ax1.set_ylabel('Depth in m')
        ax1.set_xlabel('Temperature in K')
        if count == 4:
            plt.legend(['Numerical', 'Analytical', 'Buffo'])
        else:
            pass
        display(fig)
        clear_output(wait=True)
        #plt.pause(0.1)

    def convergence_loop(self, t, T_km1, S_km1, phi_km1, Buffo=False, Stefan=False, temp_grad=None):
        # set initial values to parameters before loop
        
        [T_err, S_err, phi_err] = reset_error_for_while_loop(self.T_tol, self.S_tol, self.phi_tol) 
        if t == 1:
            [T_initial , T_km1,  T_prev,S_initial, S_prev, S_km1, phi_initial, phi_prev, phi_km1] = \
            initialize_statevariables(self.T, self.S, self.phi)
            temp_grad = self.temp_grad
        else:
            [T_initial, T_prev,S_initial,  S_prev,phi_initial,  phi_prev] = \
            set_statevariables(T_km1, S_km1, phi_km1)
        T_source = np.zeros(self.nz)
        counter = 0

        # Loop until values converge
        while T_err > self.T_tol: # HACK: or S_err > S_tol or phi_err > phi_tol:
            counter = counter + 1
            if counter > 1 :
                [T_prev, S_prev, phi_prev] = define_previous_statevariable(T_km1, S_km1, phi_km1)

            H_k = update_enthalpy(T_km1, S_km1, phi_km1, self.nz)
            H_solid = update_enthalpy_solid_state( S_km1, self.nz, liq_rel=self.liq_rel) 
            phi_k = update_liquid_fraction( T_km1,S_km1, phi_km1, H_k, H_solid, self.nz, Stefan= self.Stefan)  # BEFORE  
            advection_diffusion = AdvectionDiffusion('temperature', T_km1, T_source, T_initial, 
                                                     phi_k, phi_initial, self.w, self.dt, self.dz, 
                                                     self.nz, self.t_passed, self.S_IC, Stefan=Stefan, Buffo=Buffo, bc_neumann=self.temp_grad)  
                
            thickness, thickness_index = locate_ice_ocean_interface(phi_k, self.dz, self.nz, Stefan = self.Stefan)              
            [T_k , X_wind_T, dt_T] = advection_diffusion.unknowns_matrix()
            S_k = S_km1
            #phi_k = update_liquid_fraction( T_k,S_k, phi_km1, H_k, H_solid, self.nz, Stefan= self.Stefan)   # AFTER
            [T_km1, S_km1, phi_km1] = overwrite_statevariables(T_k, S_k, phi_k)
            
            self.thickness_index_total[t] = thickness_index
            if counter > 0:
                [T_err, T_err_full, S_err, S_err_full, phi_err, phi_err_full] = \
                    compute_error_for_convergence(T_k, T_prev, S_k, S_prev, phi_k, phi_prev)

            return (T_k, T_prev, S_k, S_prev, phi_k, phi_prev, 
                    H_k, H_solid, thickness, thickness_index, 
                    T_km1, S_km1, phi_km1)
    
    def temperature_gradient(self, phi):
        # computes the critical depth of the phase transition interface
        # for the given time step
        # returns the temperature gradient

        critical_depth = 0
        for i in range(self.nz-1):
            if phi[i] <= self.phi_c:
                critical_depth = i
                break
            else:
                pass

        self.critical_depth = self.critical_depth + (critical_depth+1)*self.dz
        return self.dz*(Tm_w - T_air_Stefan)/(self.critical_depth)

    def temperature_profile_analytical(self):
        fig1,(ax1) = plt.subplots(1,1, figsize=(10,10))
        #self.all_T, self.all_S, self.all_phi, self.all_H, self.all_H_solid, self.all_w , self.all_thick,self.all_t_passed = [[] for i in range(8)]
        T_km1, S_km1, phi_km1 = np.array([]), np.array([]), np.array([])
        T_km1_buffo, S_km1_buffo, phi_km1_buffo = np.array([]), np.array([]), np.array([])
        count = 0
        # loop over x iterations
        for t in range(1,self.iter_max):

            # loop until values converge
            if self.Buffo is True:
                T_k_buffo, T_prev_buffo, S_k_buffo, S_prev_buffo, phi_k_buffo, phi_prev_buffo, H_k_buffo, H_solid_buffo, thickness_buffo, thickness_index_buffo, T_km1_buffo, S_km1_buffo, phi_km1_buffo = \
                self.convergence_loop(t, T_km1_buffo, S_km1_buffo, phi_km1_buffo, Buffo=True, temp_grad=self.temp_grad)
            else:
                pass

            T_k, T_prev, S_k,S_prev, phi_k, phi_prev, H_k, H_solid, thickness, thickness_index, T_km1, S_km1, phi_km1 = \
            self.convergence_loop(t, T_km1, S_km1, phi_km1, Stefan=True, temp_grad=self.temp_grad)

            # computes total time passed
            self.t_passed= t_total(self.t_passed, self.dt)

            # store all results
            [self.all_T, self.all_S, self.all_phi, self.all_H, self.all_H_solid, self.all_w , self.all_thick,self.all_t_passed] = \
                store_results( T_k, S_k, phi_k, H_k, H_solid, self.w, thickness, 
                              self.t_passed, self.all_T, self.all_S, self.all_phi, 
                              self.all_H, self.all_H_solid, self.all_w, self.all_thick, 
                              self.all_t_passed, t)
            if self.bc_condition == "Neumann":
                self.temp_grad = self.temperature_gradient(phi_k)     # T_bottom - T_top / depth according to Buffo
            else:
                self.temp_grad= None
            # depth
            if self.S_IC=="S0":
                self.depth_stefan_all = stefan_problem(self.all_t_passed) 
                depth_stefan_t = self.depth_stefan_all[t]
                T_Stefan = stefan_temperature(depth_stefan_t, self.t_passed, self.dz, self.nz)
            else:
                self.depth_stefan_all = stefan_problem_twophase(self.all_t_passed) 
                depth_stefan_t = self.depth_stefan_all[t]
                T_Stefan, C_stefan = stefan_temperature_twophase(depth_stefan_t, self.t_passed, self.dz, self.nz)
            error_depth_t = (np.absolute(depth_stefan_t) - np.absolute(self.all_thick[t]))
               
            # compute temperature error and store temperature differences
            self.error_temperature[:thickness_index] = \
                np.absolute(np.absolute(T_k[:thickness_index])-np.absolute(T_Stefan[:thickness_index]))
            self.error_temperature_sum[t] = np.sum(self.error_temperature)
            self.error_temperature_sum_weighted[t] = self.error_temperature_sum[t]/thickness_index
            self.T_Stefan_diff.append(T_Stefan - self.T_Stefan_prev)
            self.T_Stefan_prev = T_Stefan
            self.T_k_diff.append(T_k - self.T_k_prev)
            self.T_k_prev = T_k
            self.T_k_list.append(T_k)
            self.T_Stefan_list.append(T_Stefan)
            if self.Buffo is True:
                self.T_k_buffo_list.append(T_k_buffo)
                self.thickness_list_Buffo.append(thickness_buffo)
            self.thickness_list.append(thickness)
            
            if t%5000 == 0:  # plot after every 500 time steps
                count = count + 1
                self.T_running(fig1, ax1, T_Stefan, T_k=T_k, T_k_buffo=T_k_buffo, count=count)

        fig1.savefig(self.folder_name + "/TemperatureProfile.png")
        self.T_k_list = np.array(self.T_k_list)
        self.T_Stefan_list = np.array(self.T_Stefan_list)
        self.T_k_diff = np.array(self.T_k_diff)
        self.T_Stefan_diff= np.array(self.T_Stefan_diff)
        self.T_k_buffo_list = np.array(self.T_k_buffo_list)

        error_norms.__init__(self, self.T_k_list, self.T_Stefan_list)

# ----------------------------------MATLAB-BUFFO-----------------------
# z_depth = 0.05
# x_axis_iter = np.arange(0,userinput.iter_max-1,1)
# T_k_ = seaicemodel_obj.T_k_list[:,int(z_depth*seaicemodel_obj.nz)]
# T_stefan_ = seaicemodel_obj.T_Stefan_list[:,int(z_depth*seaicemodel_obj.nz)]
# T_err = (np.abs(T_k_ - T_stefan_))
# index = z_depth*seaicemodel_obj.nz
# df = pd.read_csv('temp_dirichlet17.csv', sep=',')
# df_depth = np.array(df[int(z_depth*seaicemodel_obj.nz)-1: int(z_depth*seaicemodel_obj.nz)]).reshape(25000,1)

# fig1,(ax1) = plt.subplots(figsize=(10,6))
# plt.grid()
# ax1.plot(x_axis_iter*seaicemodel_obj.dt/3600, T_k_, 'k',label='Numerical Temperature')
# ax1.plot(x_axis_iter*seaicemodel_obj.dt/3600, T_stefan_, 'r--',label='Analytical Temperature')
# ax1.plot(x_axis_iter*seaicemodel_obj.dt/3600,df_depth[:24999], label='Buffo' )
# ax1.set_xlabel("t in hours")
# ax1.set_ylabel("Temperature in K")
# #ax1.legend()
# ax1.set_title("Temperature evolution at {}m".format(z_depth))
# color = 'gray'
# ax3 = ax1.twinx()
# ax3.plot(x_axis_iter*seaicemodel_obj.dt/3600, T_err, color=color, label='|$T_{Numerical}$-$T_{Analytical}$|',linestyle='dashed', linewidth=1)
# ax3.tick_params(axis='y', labelcolor=color)
# ax1.legend()
# fig1.tight_layout()
#-----------------------------------Interpolation------------------------------------------------------
# #temp_01 = pd.read_csv("Temperature_0.01\Temperature0.01.csv", sep=',')
# temp_001 = pd.read_csv("Temperature_0.001\Temperature0.001.csv", sep=',')
# temp_005 = pd.read_csv("Temperature_0.005\Temperature0.005.csv", sep=',')
# i = 20000
# temp_001  = temp_001.drop(['Unnamed: 0'], axis=1)
# temp_005 = temp_005.drop(['Unnamed: 0'], axis=1)
# #temp_01_i100 = np.array(temp_01[i:i+1])
# temp_005_i100 = np.array(temp_005[i:i+1])
# temp_001_i100 = np.array(temp_001[i:i+1])
# x_new = [x for x in range(temp_001_i100.shape[1])]
# #y_001 = np.interp(x_new, [x for x in range(temp_001_i100.shape[1])], temp_001_i100[0])
# y_005 = np.interp(x_new,[x for x in range(temp_005_i100.shape[1])], temp_005_i100[0])

# plt.grid()
# plt.plot(x_new, y_005, label='dz=0.001')
# plt.plot(x_new, temp_001_i100[0], label='dz=0.005')
# plt.xlabel('No. of nodes')
# plt.ylabel('Temperature')
# plt.legend()
# plt.title('Temperature at iter=20000')
# plt.show()

# plt.grid()
# plt.plot([x for x in range(temp_005_i100.shape[1])], temp_005_i100[0], label='dz=0.005')
# plt.plot(x_new, temp_001_i100[0], label='dz=0.001')
# plt.xlabel('No. of nodes')
# plt.ylabel('Temperature')
# plt.legend()
# plt.title('Temperature at iter=20000, wo_interp')
# plt.show()
#-------------------------------ErrorNorm------------------------------------------------------
# from calculate_error import error_norms
# err_dz = error_norms(y_005,temp_001_i100[0])
# diff_005_001 = err_dz.numerical_analytical_diff()
# one_norm_005_001 = err_dz.two_norm(diff_005_001)
# inf_norm_005_001 = err_dz.infinity_norm(diff_005_001)

# plt.grid()
# plt.plot(x_new, one_norm_005_001, label='one norm')
# plt.plot(x_new, inf_norm_005_001, label='two norm')
# plt.xlabel('No. of nodes')
# plt.ylabel('|Temp_005 - Temp_001|')
# plt.legend()
# plt.title('Temperature at iter=20000, wo_interp')
# plt.show()

#-----------------------------CompareSolver----------------------------------------------------
# def TDMA_V3(a,b,c,d):
#     n = len(d)
#     w= np.zeros(n-1,float)
#     g= np.zeros(n, float)
#     p = np.zeros(n,float)
    
#     w[0] = c[0]/b[0]
#     g[0] = d[0]/b[0]

#     for i in range(1,n-1):
#         w[i] = c[i]/(b[i] - a[i-1]*w[i-1])
#     for i in range(1,n):
#         g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1])
#     p[n-1] = g[n-1]
#     for i in range(n-1,0,-1):
#         p[i-1] = g[i-1] - w[i-1]*p[i]
#     return p

# lower = np.array([1,2,3])
# upper = np.array([2,1,1])
# main = np.array([4,7,6,5])
# rhs = np.array([1,1,1,1])
# solved = TDMA_V3(a=lower,b=main,c=upper,d=rhs)

# def solver_compare(lower, main, upper):
#     rhs = np.array([1,1,1,1])
#     diff_crit = np.abs(main[:-1]) - (np.abs(lower) + np.abs(upper)) 
#     solved = TDMA_V3(a=lower,b=main,c=upper,d=rhs)
#     fin_mat = np.diag(lower, k=1) + np.diag(upper, k=-1) + np.diag(main, k=0)
#     solved_numpy = np.linalg.solve(fin_mat,rhs)
#     solver_diff = solved_numpy - solved
#     return diff, solver_diff      
