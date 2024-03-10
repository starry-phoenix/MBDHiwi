import numpy as np
import os
import sys
import datetime
import hickle as hkl


def set_up_matrices(iter_max, nz):
    """
    Initialize matrices with one column for each time step
    """
    all_T = np.zeros([iter_max, nz])
    all_S_sw = np.zeros([iter_max, nz])
    all_phi = np.zeros([iter_max, nz])
    all_H = np.zeros([iter_max, nz])
    all_H_solid = np.zeros([iter_max, nz])
    all_w = np.zeros([iter_max, nz])
    all_thick = np.zeros([iter_max])
    all_t_passed = np.zeros(iter_max)
    return all_T, all_S_sw, all_phi, all_H, all_H_solid, all_w, all_thick, all_t_passed


def store_results(
    T,
    S_sw,
    phi,
    H,
    H_solid,
    w,
    thickness,
    t_passed,
    all_T,
    all_S_sw,
    all_phi,
    all_H,
    all_H_solid,
    all_w,
    all_thick,
    all_t_passed,
    t
):
    all_T[t, :] = T
    all_S_sw[t, :] = S_sw
    all_phi[t, :] = phi
    all_H[t, :] = H
    all_H_solid[t, :] = H_solid
    all_w[t, :] = w
    all_thick[t] = thickness
    all_t_passed[t] = t_passed
    return all_T, all_S_sw, all_phi, all_H, all_H_solid, all_w, all_thick, all_t_passed


def create_directory():
    # Parent Directory path
    parent_dir = "/Users/asimson/work/GIT_projects/seaicemodel/data"  # "C:\\Users\\Anna Lara\\Documents\\05_GIT\\Seaicemodel\\Data\\"
    # Directory

    directory = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    # Path
    path = os.path.join(parent_dir, directory)
    os.mkdir(path)
    print("Directory '% s' created" % directory)
    sys.path.append(path)


def save_results(
    all_T, all_S, all_phi, all_H, all_H_solid, all_w, all_thick, all_t_passed
):
    """
    save results as hkl files
    """
    # go to path
    os.chdir(sys.path[-1])
    # save data
    hkl.dump(all_T, "temperature.hkl", compression="gzip")
    print("saved temperature")
    hkl.dump(all_phi, "liquid_fraction.hkl", compression="gzip")
    print("saved liquid fraction")
    hkl.dump(all_S, "salinity.hkl", compression="gzip")
    print("saved salinity")
    hkl.dump(all_H, "enthalpy.hkl", compression="gzip")
    print("saved enthalpy")
    hkl.dump(all_H_solid, "enthalpy_solid.hkl", compression="gzip")
    print("saved solid state enthalpy")
    hkl.dump(all_w, "velocity.hkl", compression="gzip")
    print("saved velocity")
    hkl.dump(all_thick, "thickness.hkl", compression="gzip")
    print("saved thickness")
    hkl.dump(all_t_passed, "time_passed.hkl", compression="gzip")
    print("saved time passed")
