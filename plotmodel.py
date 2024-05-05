import matplotlib.pyplot as plt
from calculate_error import error_norms
from seaicemodel import SeaIceModel
import numpy as np
import userinput
import time
import pandas as pd
from model.constants_real import Tm_w, S_sw
np.seterr(divide='ignore', invalid='ignore')

class PlotModel(SeaIceModel):
    def __init__(self) -> None:
        SeaIceModel.__init__(self)
        
    def error_analytical_numerical(self):
        # Calculates the errors of Numerical and Analytical using error norms one, two and infinity
        # norm(|T_k_list - T_Stefan_list|)
        
        num_ana_temperature_diff = self.numerical_analytical_diff()
        self.T_k_Stefan_diff_L1norm = self.one_norm(num_ana_temperature_diff)
        self.T_k_Stefan_diff_infnorm = self.infinity_norm(num_ana_temperature_diff)
        self.T_k_Stefan_diff_L2norm = self.two_norm(num_ana_temperature_diff)

        # Calculates the errors of Analytical temperature differences between two consecutive iterations
        # norm(|T_Stefan[i+1] - T_Stefan[i]|)
        self.T_Stefan_diff_infnorm = self.infinity_norm(self.T_Stefan_diff)   
        self.T_Stefan_diff_L2norm = self.two_norm(self.T_Stefan_diff)
        self.T_Stefan_diff_L1norm = self.one_norm(self.T_Stefan_diff)

        # Calculates the errors of Numerical temperature differences between two consecutive iterations
        # norm(|T_k[i+1] - T_k[i]|)
        self.T_k_diff_infnorm = self.infinity_norm(self.T_k_diff)
        self.T_k_diff_L2norm = self.two_norm(self.T_k_diff)
        self.T_k_diff_L1norm = self.one_norm(self.T_k_diff)
    
    def phi_slope(self,iteration):
        phi = self.phi_k_list[iteration]
        # T = self.T_k_list[iteration]
        # phi_diff = np.roll(phi,-1)[:-1]-phi[:-1] 
        # T_diff = np.roll(T,-1)[:-1]- T[:-1] 
        # slope = phi_diff/T_diff
        # slope = np.round(np.nan_to_num(slope, posinf=0.0), 1)
        mush_idx = np.intersect1d(np.where(phi<1.0), np.where(phi>0.0))
        if mush_idx.size == 0:
            mush_idx = np.where(phi==1.0)[0]
        return mush_idx

    def plot_error_temp_diff(self, zoom_x, savefig="True"):
        x_axis_iter = np.arange(0,zoom_x,1)
        fig, ax1 = plt.subplots(figsize=(10,6))
        plt.grid()
        #ax1.plot(x_axis_iter*self.dt/3600, self.T_k_Stefan_diff_infnorm[x_axis_iter], 'k--',label='$L_\infty$')
        ax1.plot(x_axis_iter*self.dt/3600, self.T_k_Stefan_diff_L2norm[x_axis_iter], 'k',label='$L_2$')
        #ax1.plot(x_axis_iter, self.T_k_Stefan_diff_L1norm[x_axis_iter], label='$Norm_1$')
        ax1.set_yscale('log')
        ax1.set_xlabel("t in hours")
        ax1.set_ylabel("$\Delta T$")
        ax1.set_title("|T_Stefan - T_k|")
        ax1.tick_params(axis='y')
        ax1.legend()
        ax1.set_yscale('log')
        ax2 = ax1.twinx() 
        color = 'teal'
        ax2.plot(x_axis_iter*self.dt/3600, self.thickness_list[:len(x_axis_iter)], label='Thickness in m', color=color, alpha=0.7, linestyle='dashed')
        ax2.set_ylabel("Depth in m", color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        #plt.xscale('log')
        #ax2.legend()
        fig.tight_layout()      
        plt.show()
        if savefig:
            fig.savefig(self.folder_name + "/Temperature_error_diff_num_ana.png")
        else:
            pass

    def plot_error_temp(self, zoom_x, norm="inf", savefig="True"):
        fig1,(ax1) = plt.subplots(figsize=(10,6))
        plt.grid()
        x_axis_iter = np.arange(0,zoom_x,1)
        #plt.figure(figsize=(10,10))
        if norm=="inf":
            ax1.plot(x_axis_iter*self.dt/3600, self.T_Stefan_diff_infnorm[x_axis_iter], 'r--',label='$L_\infty$ Analytical')
            ax1.plot(x_axis_iter*self.dt/3600, self.T_k_diff_infnorm[x_axis_iter], 'k',label='$L_\infty$ Numerical')
        elif norm=="2":
            ax1.plot(x_axis_iter*self.dt/3600, self.T_Stefan_diff_L2norm[x_axis_iter], 'r--',label='$L_2$ Analytical')
            ax1.plot(x_axis_iter*self.dt/3600, self.T_k_diff_L2norm[x_axis_iter], 'k',label='$L_2$ Numerical')
        else:
            ax1.plot(x_axis_iter*self.dt/3600, self.T_Stefan_diff_L1norm[x_axis_iter], 'r--',label='$L_1$ Analytical')
            ax1.plot(x_axis_iter*self.dt/3600, self.T_k_diff_L1norm[x_axis_iter], 'k',label='$L_1$ Numerical')
        ax1.set_xlabel("t in hours")
        ax1.set_ylabel("$\Delta T$")
        ax1.set_title("Temperature Error of Analytical and Numerical Results")
        ax1.set_yscale('log')
        ax1.legend()
        color = 'teal'
        ax3 = ax1.twinx()
        ax3.plot(x_axis_iter*self.dt/3600, self.thickness_list[:len(x_axis_iter)], color=color, alpha=0.7)
        #ax3.set_ylabel("Thickness in m")
        ax3.tick_params(axis='y', labelcolor=color)

        fig1.tight_layout()
        plt.show()
        if savefig:
            fig1.savefig(self.folder_name + "/Temperature_error_num&ana.png")
        else:
            pass
    
    def plot_depth_over_time(self):
        x_axis_iter = np.arange(0,userinput.iter_max-1,1)
        depth_mush = np.append(np.arange(0,1,0.01),1.0)
        mush_list_y1 = np.array([[depth_mush[self.phi_slope(i)[0]], depth_mush[self.phi_slope(i)[-1]]] for i in x_axis_iter])
        #mush_list_y2 = [depth_mush[self.phi_slope(i)[-1]] for i in x_axis_iter]
        depth = np.array(self.thickness_list_Buffo[:len(x_axis_iter)])
        plt.figure(figsize=(10,6))
        plt.grid()
        plt.plot(x_axis_iter*self.dt/3600, self.thickness_list[:len(x_axis_iter)], 'r--',label='Numerical Depth')
        plt.plot(x_axis_iter*self.dt/3600, self.depth_stefan_all[:len(x_axis_iter)], 'k',label='Analytical Depth')
        if self.Buffo is True:
            plt.plot(x_axis_iter*self.dt/3600, self.thickness_list_Buffo[:len(x_axis_iter)], 'b-.',alpha=0.5,label='Buffo Depth')
        plt.fill_between(x_axis_iter*self.dt/3600, mush_list_y1[:,0], mush_list_y1[:,1], color='gray', alpha=0.2, label='Mushy Layer')
        plt.xlabel("t in hours")
        plt.ylabel("Depth in m")
        plt.legend()
        plt.title("Numerical Depth Vs Analytical Depth")
        plt.savefig(self.folder_name +"/Numerical_Analytical_Depth.png")
        plt.show()

    def plot_temperature(self, z_depth, savefig=True, Buffo_matlab=False):
        x_axis_iter = np.arange(0,userinput.iter_max-1,1)
        #x_axis_iter = np.arange(0,22970,1)
        T_k_ = self.T_k_list[:,int(z_depth*self.nz)]
        T_stefan_ = self.T_Stefan_list[:,int(z_depth*self.nz)]
        T_err = (np.abs(T_k_ - T_stefan_))
        if self.Buffo is True:
            T_k_buffo_ = self.T_k_buffo_list[:,int(z_depth*self.nz)]
            T_err_buffo = (np.abs(T_k_ - T_k_buffo_))
        index = z_depth*self.nz
        if Buffo_matlab is True:
            df = pd.read_csv('MatlabData/temp_dirichletSalinity01.csv', sep=',')
            df_depth = np.array(df[int(z_depth*self.nz)-1: int(z_depth*self.nz)]).reshape(25000,1)
        
        fig1,(ax1) = plt.subplots(figsize=(10,6))
        plt.grid()
        ax1.plot(x_axis_iter*self.dt/3600, T_k_[x_axis_iter], 'r--',label='Numerical Temperature')
        ax1.plot(x_axis_iter*self.dt/3600, T_stefan_[x_axis_iter], 'k',label='Analytical Temperature')
        if self.Buffo is True:
            ax1.plot(x_axis_iter*self.dt/3600, T_k_buffo_[x_axis_iter], 'b--',alpha=0.5,label='Buffo')
        if Buffo_matlab is True:
            ax1.plot(x_axis_iter*self.dt/3600,df_depth[:24999], 'b',alpha=0.2,label='Buffo-matlab' )
        ax1.set_xlabel("t in hours")
        ax1.set_ylabel("Temperature in K")
        #ax1.legend()
        ax1.set_title("Temperature evolution at {}m".format(z_depth))
        color = 'gray'
        ax3 = ax1.twinx()
        ax3.plot(x_axis_iter*self.dt/3600, T_err[x_axis_iter], color=color, label='|$T_{Numerical}$-$T_{Analytical}$|',linestyle='dashed', linewidth=1)
        ax3.tick_params(axis='y', labelcolor=color)
        ax3.set_ylabel('|T_k - T_Stefan|', color=color)
        ax1.legend()
        fig1.tight_layout()
        if savefig:
            fig1.savefig(self.folder_name + "/Temperature evolution at"+ str(z_depth) +"m.png")
        else:
            pass
        plt.show()

        df_temp = pd.DataFrame(self.T_k_list)
        df_temp.to_csv(self.folder_name + "/Temperature" + str(self.dz) + '.csv')

    def plot_phi(self, timestep, savefig=True):
        x_axis_iter = np.arange(0,101,1)
        #x_axis_iter = np.arange(0,22970,1)
        index = timestep*3600/self.dt
        phi_k = self.phi_k_list[int(index)]
        if self.Buffo is True:
            phi_buffo = self.phi_buffo_list[int(index)]
        
        fig1,(ax1) = plt.subplots(figsize=(10,6))
        plt.grid()
        ax1.plot(x_axis_iter*self.dt/3600, phi_k, 'r--',label='Numerical Temperature')
        if self.Buffo is True:
            ax1.plot(x_axis_iter*self.dt/3600, phi_buffo, 'b--',alpha=0.5,label='Buffo')
        ax1.set_xlabel("depth in nodes")
        ax1.set_ylabel("Phi")
        #ax1.legend()
        ax1.set_title("Temperature evolution at {}H".format(timestep))
        color = 'gray'
        ax1.legend()
        fig1.tight_layout()
        if savefig:
            fig1.savefig(self.folder_name + "/Liquid Fraction evolution at"+ str(timestep) +"s.png")
        else:
            pass
        plt.show()

    def plot_salinity(self, z_depth, savefig=True):
        x_axis_iter = np.arange(0,userinput.iter_max-1,1)
        index = z_depth*self.nz
        S_k = self.S_k_list[:,int(index)]
        if self.Buffo is True:
            S_buffo = self.S_buffo_list[:,int(index)]
        fig1,(ax1) = plt.subplots(figsize=(10,6))
        plt.grid()
        ax1.plot(x_axis_iter*self.dt/3600, S_k[x_axis_iter], 'r--',label='Numerical Temperature')
        if self.Buffo is True:
            ax1.plot(x_axis_iter*self.dt/3600, S_buffo[x_axis_iter], 'b--',alpha=0.5,label='Buffo')
        ax1.set_xlabel("t in hours")
        ax1.set_ylabel("Salinity in ppt")
        #ax1.legend()
        ax1.set_title("Salinity evolution at {}m".format(z_depth))
        color = 'gray'
        ax1.legend()
        fig1.tight_layout()
        if savefig:
            fig1.savefig(self.folder_name + "/Salinity evolution at"+ str(z_depth) +"m.png")
        else:
            pass
        plt.show()

    def plot_enthalpy(self, timestep, savefig=False):
        x_axis_iter = np.arange(0,userinput.iter_max-1,1)
        index = timestep*3600/self.dt
        mush = self.phi_slope(int(index))
        phi = self.phi_k_list[int(index)]
        H_k = self.H_k_list[int(index)]
        H_solid = self.H_solid_list[int(index)]
        T_k = self.T_k_list[int(index)]
        H = (H_k - H_solid)/334774
        T_melt = Tm_w 
        T_interface = self.Temp_interface[int(index)]
        fig1,(ax1) = plt.subplots(figsize=(10,6))
        plt.grid()
        ax1.plot(T_k,phi, 'r--',label='Phi')
        ax1.fill_betweenx(H, T_k[mush][0], T_k[mush][-1], color='gray', alpha=0.2, label='Mushy Layer')
        ax1.set_xlabel("Temperature in K")
        ax1.set_ylabel("Liquid Fraction $\phi$")
        ax1.set_title("Liquid Fraction Vs Temperature at {}h".format(timestep))
        color = 'gray'
        ax3 = ax1.twinx()
        ax3.axvline(T_melt, color='b', linestyle='dashed', label='T_melt:'+str(round(T_melt,2)))
        ax3.axvline(T_interface, color='g', linestyle='dashed', label='T_interface:' + str(round(T_interface,2)))
        fig1.tight_layout()
        ax1.legend(loc=0)
        ax3.legend(loc=1)
        if savefig:
            fig1.savefig(self.folder_name + "/LiqVsTemp evolution at"+ str(timestep) +"h.png")
        else:
            pass
        plt.show()


#from model.constants_real import k_i, k_w

# dz, dt, iteration = plot_seaicemodel_obj.dz , plot_seaicemodel_obj.dt, plot_seaicemodel_obj.iter_max
# T = plot_seaicemodel_obj.T_k_list
# k_x = lambda phi_x: phi_x*k_w + (1-phi_x)*k_i
# alpha = dt/dz
# Q_diff = np.zeros((iteration, 99))
# H_diff = np.zeros((iteration-1, 100))
# H = plot_seaicemodel_obj.H_k_list

# def flux_i_n(i,n):
#     R_i = 0.5*(dz/k_x(phi[n,i]) + dz/k_x(phi[n,i+1]))
#     Q_i = (T[n,i] - T[n,i+1])/R_i
    
#     return Q_i

# for t in range(iteration-1):
#     Q_diff[t] = [(flux_i_n(i,t) - flux_i_n(i-1,t))/dz for i in range(1,100)]

# H_diff = (H[1:,1:] - H[:-1,1:])/dt  

# import matplotlib.pyplot as plt 

# ax_time = np.arange(1,24998,1)
# ax_size = np.arange(1,100)
# t_ = 10
# plt.plot(ax_size, H_diff[t_,:-1], label='Enthalpy diff')
# plt.plot(ax_size, Q_diff[t_], label='Flux diff')
# plt.xlabel('space dz nodes')
# plt.legend()
# plt.show()
