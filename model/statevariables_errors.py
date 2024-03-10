import numpy as np

# Function to set state variables
def set_statevariables(T_km1, S_km1, phi_km1):
    T_prev, T_initial = T_km1, T_km1  # Set previous T as current T for initialization, # Store initial T for later use
    S_prev, S_initial = S_km1, S_km1  # Set previous S as current S for initialization, # Store initial S for later use
    phi_prev, phi_initial = phi_km1, phi_km1  # Set previous phi as current phi for initialization, # Store initial phi for later use

    return T_initial, T_prev, S_initial, S_prev, phi_initial, phi_prev


# Function to define previous state variables
def define_previous_statevariable(T_km1, S_km1, phi_km1):
    T_prev = T_km1  # Set previous T as current T for iteration
    S_prev = S_km1  # Set previous S as current S for iteration
    phi_prev = phi_km1  # Set previous phi as current phi for iteration
    return T_prev, S_prev, phi_prev


# Function to overwrite state variables
def overwrite_statevariables(T_k, S_k, phi_k):
    T_km1 = T_k  # Set current T as previous T for next iteration
    S_km1 = S_k  # Set current S as previous S for next iteration
    phi_km1 = phi_k  # Set current phi as previous phi for next iteration
    return T_km1, S_km1, phi_km1


# Function to compute errors for convergence
def compute_error_for_convergence(T_k, T_prev, S_k, S_prev, phi_k, phi_prev):
    T_err = np.max(abs(
        (T_k[1:-1] - T_prev[1:-1]))
    )  # Compute maximum T error for convergence check
    T_err_full = abs(T_k - T_prev)  # Compute full T error for convergence check
    S_err = np.max(
        (S_k[1:-1] - S_prev[1:-1])
    )  # Compute maximum S error for convergence check
    S_err_full = S_k - S_prev  # Compute full S error for convergence check
    phi_err = np.max(
        (phi_k[1:-1] - phi_prev[1:-1])
    )  # Compute maximum phi error for convergence check
    phi_err_full = phi_k - phi_prev  # Compute full phi error for convergence check
    return T_err, T_err_full, S_err, S_err_full, phi_err, phi_err_full


# Function to reset errors for while loop
def reset_error_for_while_loop(T_tol, S_tol, phi_tol):
    T_err = 1 + T_tol  # Set initial T error to value greater than tolerance
    S_err = 1 + S_tol  # Set initial S error to value greater than tolerance
    phi_err = 1 + phi_tol  # Set initial phi error to value greater than tolerance
    return T_err, S_err, phi_err


# Function to initialize state variables
def initialize_statevariables(T, S, phi):
    phi_initial = phi  # Store initial phi for later use
    phi_km1 = phi  # Set previous phi as current phi for iteration
    phi_prev = phi  # Set previous phi as current phi for initialization
    T_initial = T  # Store initial T for later use
    T_prev = T  # Set previous T as current T for initialization
    T_km1 = T  # Set current T as previous T for iteration
    S_initial = S  # Store initial S for later use
    S_km1 = S  # Set current S as previous S for initialization
    S_prev = S
    return (
        T_initial,
        T_km1,
        T_prev,
        S_initial,
        S_km1,
        S_prev,
        phi_initial,
        phi_km1,
        phi_prev,
    )
