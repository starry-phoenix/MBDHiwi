import matplotlib.pyplot as plt
from calculate_error import error_norms
from seaicemodel import SeaIceModel
import numpy as np
import userinput
import time
import pandas as pd

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
        plt.figure(figsize=(10,6))
        plt.grid()
        plt.plot(x_axis_iter*self.dt/3600, self.thickness_list[:len(x_axis_iter)], 'r--',label='Numerical Depth')
        plt.plot(x_axis_iter*self.dt/3600, self.depth_stefan_all[:len(x_axis_iter)], 'k',label='Analytical Depth')
        if self.Buffo is True:
            plt.plot(x_axis_iter*self.dt/3600, self.thickness_list_Buffo[:len(x_axis_iter)], 'b-.',alpha=0.5,label='Buffo Depth')
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
