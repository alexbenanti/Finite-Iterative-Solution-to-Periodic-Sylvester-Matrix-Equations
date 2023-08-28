import numpy as np 
import matplotlib.pyplot as plt 

def Q_calculator_backward(A,B,C,X,periodicity):
    Q = []
    for j in range(periodicity):
        Q_j = C[j] - A[j]@X[(j+1)%3] - X[j]@B[j]
        Q.append(Q_j)
    return Q

def R_calculator_backward(A,Q,B,periodicity):
    R = []
    R_norm = 0
    for j in range(periodicity):
        R_j = A[j-1].T@Q[j-1] + Q[j]@B[j].T
        R_norm += np.linalg.norm(R_j)**2
        R.append(R_j)
    return R, R_norm

def initialize_P(R):
    return -np.array(R)

def alpha_calculator(A, B, R, P, periodicity):
    trace = 0
    norm = 0
    for j in range(periodicity):
        trace+= np.trace(P[j].T@R[j])
        norm+= np.linalg.norm(A[j]@P[(j+1)%3] + P[j]@B[j])**2
    return trace/norm

def update_X(X_current, P, alpha):
    return X_current + alpha*P

def P_updator(P, R, new_norm, old_norm, peridocity):
    return -np.array(R) + (new_norm/old_norm)*P

def error_calculator(X_prime, X_true, periodicity):
    return sum([np.linalg.norm(X_prime[j]-X_true[j])**2 for j in range(periodicity)])


def plot_error_b(deltas):
    iterations = np.linspace(0,len(deltas), len(deltas))
    zeros = np.zeros(len(iterations))
    plt.plot(iterations,deltas,'-b')
    plt.plot(iterations,zeros, '-g')
    plt.title('Error vs. Number of Iterations (Forward Equations)')
    plt.xlabel('Iteration Number ($k$)')
    plt.ylabel('Error')
    plt.show()