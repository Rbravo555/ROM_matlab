# Discrete Empirical Interpolation Method
# Based on paper from Chaturantabut & Sorensen, 2010

import random
import numpy as np
from matplotlib import pyplot as plt
from randomized_singular_value_decomposition import RandomizedSingularValueDecomposition


def L2_difference(Original, Approximation):
    percentual_diff = (np.linalg.norm( Approximation -  Original ))   /  (np.linalg.norm(Original)) * 100
    return percentual_diff


def DEIM(Basis):
    #find first point
    U = Basis[:,0].reshape(Basis.shape[0], -1)
    z = np.zeros(U.shape)
    P = z.copy()
    indexes = np.argmax( np.abs(Basis[:,0]) )
    P[indexes] = 1

    #find next points
    for i in range(1,Basis.shape[1]):
        c = np.linalg.solve(P.T @ U , P.T @ Basis[:,i] )
        residual = Basis[:,i] - U @ c
        U = np.c_[U,Basis[:,i]]
        index_i = np.argmax( np.abs(residual) )
        indexes = np.r_[indexes, index_i]
        P = np.c_[P, z]; P[index_i,i] = 1

    return indexes


def compare_random_gappy_vs_DEIM():

    S = np.load('SnapshotsMatrix.npy')
    u,s,v,error = RandomizedSingularValueDecomposition().Calculate(S, 1e-4)
    number_of_modes = int(len(s))
    print('number of modes taken:', number_of_modes)

    #randomly selecting the same number of points as the number of columns in the basis
    gappy_points = random.sample(range(0, np.shape(S)[0]), number_of_modes)
    print('gappy_points:', gappy_points)
    u_gappy = u[gappy_points, :]
    S_gappy = S[gappy_points, :]
    S_reconstructed_gappy = u @ (np.linalg.pinv(u_gappy) @ S_gappy)

    #using DEIM to select the best points
    DEIM_points = DEIM(u)
    print('DEIM_points:', DEIM_points)
    u_DEIM = u[DEIM_points, :]
    S_DEIM = S[DEIM_points, :]
    S_reconstructed_DEIM = u @ (np.linalg.pinv(u_DEIM) @ S_DEIM)


    print('Difference in the reconstruction using ', len(s),'random points', L2_difference(S, S_reconstructed_gappy)  )
    print('Difference in the reconstruction using ', len(s),'DEIM points', L2_difference(S, S_reconstructed_DEIM)  )


if __name__=='__main__':
    compare_random_gappy_vs_DEIM()

