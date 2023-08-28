import numpy as np
from backward_utils import *

def backward_algorithm(A,B,C,X_initial, X_true, periodicity, tolerance, max_iter):
    Q_init = Q_calculator_backward(A,B,C,X_initial, periodicity)
    true_norms_sum = sum([np.linalg.norm(X_true[j])**2 for j in range(periodicity)])
    initial_error = error_calculator(X_initial, X_true, periodicity)
    deltas = [np.sqrt(initial_error/true_norms_sum)]
    R, R_norm_last = R_calculator_backward(A,Q_init, B, periodicity)
    P = initialize_P(R)
    X_current  = X_initial
    for i in range(max_iter):
        if all(np.linalg.norm(R_j)<tolerance for R_j in R):
             break
        alpha = alpha_calculator(A,B,R,P, periodicity)
        X_new = update_X(X_current, P, alpha)
        Q_new = Q_calculator_backward(A,B,C,X_new, periodicity)
        R_new, R_new_norm  = R_calculator_backward(A,Q_new, B, periodicity)
        P_new = P_updator(P, R_new, R_new_norm, R_norm_last, periodicity)
        current_error = error_calculator(X_new, X_true, periodicity)
        deltas.append(np.sqrt(current_error/true_norms_sum))
        P = P_new
        R = R_new
        R_norm_last = R_new_norm
        X_current = X_new
    
    return deltas