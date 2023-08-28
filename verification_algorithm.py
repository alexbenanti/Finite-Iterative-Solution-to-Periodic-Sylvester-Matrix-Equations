import numpy as np
from forward_utils import *

def iterative_algorithm(A,B,C,X_initial, X_true, periodicity, tolerance, max_iter):
    Q_init = Q_calculator(A,B,C,X_initial, periodicity)
    R, R_norm_last = R_calculator(A,Q_init, B, periodicity)
    P = initialize_P(R)
    X_current  = X_initial
    X0_00 = [X_current[0][0,0]]
    X0_01 = [X_current[0][0,1]]
    X0_10 = [X_current[0][1,0]]
    X0_11 = [X_current[0][1,1]]
    X1_00 = [X_current[1][0,0]]
    X1_01 = [X_current[1][0,1]]
    X1_10 = [X_current[1][1,0]]
    X1_11 = [X_current[1][1,1]]
    X2_00 = [X_current[2][0,0]]
    X2_01 = [X_current[2][0,1]]
    X2_10 = [X_current[2][1,0]]
    X2_11 = [X_current[2][1,1]]
    for i in range(max_iter):
        if all(np.linalg.norm(R_j)<tolerance for R_j in R):
             break
        alpha = alpha_calculator(A,B,R,P, periodicity)
        X_new = update_X(X_current, P, alpha)
        Q_new = Q_calculator(A,B,C,X_new, periodicity)
        R_new, R_new_norm  = R_calculator(A,Q_new, B, periodicity)
        P_new = P_updator(P, R_new, R_new_norm, R_norm_last, periodicity)
        P = P_new
        R = R_new
        R_norm_last = R_new_norm
        X0_00.append(X_new[0][0,0])
        X0_01.append(X_new[0][0,1])
        X0_10.append(X_new[0][1,0])
        X0_11.append(X_new[0][1,1])
        X1_00.append(X_new[1][0,0])
        X1_01.append(X_new[1][0,1])
        X1_10.append(X_new[1][1,0])
        X1_11.append(X_new[1][1,1])
        X2_00.append(X_new[2][0,0])
        X2_01.append(X_new[2][0,1])
        X2_10.append(X_new[2][1,0])
        X2_11.append(X_new[2][1,1])
        X_current = X_new
    
    return X0_00,X0_01,X0_10,X0_11,X1_00,X1_01,X1_10,X1_11,X2_00,X2_01,X2_10,X2_11, X_current
        