"""
Set up time parameters
"""
import numpy as np
import userinput
### Initialize time step, starting time, end time

### Computes total time passed


def t_total(t_passed, dt):
    """
    computes total time passed
    based on current time step dt and, total time of previous time step
    """
    return t_passed + dt


### Computes maximum number of iterations based on start and end time and time step


def set_up_iter(iter_max):
    """
    Set up maximum number of iterations (not related to time)
    Set up first time step
    """

    dt = userinput.dt
    t_passed0 = 0
    iter_max = iter_max
    return iter_max, dt, t_passed0


# def CFL_time_step(w,dt,dz,nz):
#     '''
#     time step advection
#     '''
#     dt = np.zeros(nz)
#     for i in range(nz):
#         dt[i] = abs(dz/w[i])
#     dt_non_zero = dt[dt !=0]
#     dt_CFL = np.amin(dt_non_zero)
#     return dt_CFL

