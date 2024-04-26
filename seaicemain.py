import time 
from plotmodel import PlotModel
# make changed in userinput.py before running this script
if __name__ == "__main__":
    start = time.time()
    plot_seaicemodel_obj = PlotModel()
    plot_seaicemodel_obj.temperature_profile_analytical()
    end = time.time()
    run_time = end-start
    plot_seaicemodel_obj.error_analytical_numerical()
    plot_seaicemodel_obj.plot_error_temp(5000, norm='inf', savefig=True)     # norm="inf" OR norm="2" OR norm="1"
    plot_seaicemodel_obj.plot_error_temp_diff(5000, savefig=True)
    plot_seaicemodel_obj.plot_depth_over_time()
    plot_seaicemodel_obj.plot_temperature(z_depth=0.1, savefig=True, Buffo_matlab=True)
    plot_seaicemodel_obj.plot_phi(timestep=10, savefig=True)
    plot_seaicemodel_obj.plot_enthalpy(timestep=10, savefig=False)