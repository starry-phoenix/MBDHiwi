import numpy as np
import sys
class error_norms:
    def __init__(self, numerical_values, analytical_values) -> None:
        self.numerical_values = numerical_values
        self.analytical_values = analytical_values
        self.iteration_count = 0
        
    def numerical_analytical_diff(self):
        if len(self.analytical_values.shape) ==1:
            self.iteration_count, self.array_len = 1, self.analytical_values.shape[0]
            self.numerical_analytical_diff = abs(self.analytical_values - self.numerical_values)
        else:
            [self.iteration_count, self.array_len] = self.analytical_values.shape
            self.numerical_analytical_diff = np.array([abs(self.analytical_values[i] - self.numerical_values[i]) for i in range(self.iteration_count)])
        return self.numerical_analytical_diff

    def one_norm(self, numerical_analytical_diff):
        if len(self.analytical_values.shape) ==1:
            one_norm_result = np.sum(numerical_analytical_diff)/self.array_len
        else:
            one_norm_result = np.array([np.sum(numerical_analytical_diff[i])/self.array_len for i in range(self.iteration_count)])
        return one_norm_result
    
    def two_norm(self, numerical_analytical_diff):
        if len(self.analytical_values.shape) ==1:
            two_norm_result = np.sqrt(np.sum(numerical_analytical_diff**2))/self.array_len
        else:
            two_norm_result = np.array([np.sqrt(np.sum(numerical_analytical_diff[i]**2))/self.array_len for i in range(self.iteration_count)])
        return two_norm_result
    
    def infinity_norm(self, numerical_analytical_diff):
        if len(self.analytical_values.shape) ==1:
            inf_norm_result = max(numerical_analytical_diff)
        else:
            inf_norm_result = np.array([max(numerical_analytical_diff[i]) for i in range(self.iteration_count)])
        return inf_norm_result