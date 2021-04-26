from SetRefereceElement import setup_reference_element
from SetDomain import setup_domain

import numpy as np
from matplotlib import pyplot as plt
from celluloid import Camera


class FOM_simulation(object):

    def __init__( self, number_of_elements=10, number_of_steps=10,element_degree=1):
        self.domain = setup_domain(number_of_elements)
        self.reference_element = setup_reference_element(element_degree)
        self.number_of_steps = number_of_steps

    def Run(self):
        self.ComputeSystemMatrix()
        self.Solve()

    def ComputeSystemMatrix(self):
        self.K = np.zeros((self.domain.number_of_nodes, self.domain.number_of_nodes))
        #element loop
        for ith_element in range(self.domain.number_of_elements):
            element_connectivity = self.domain.connectivity[ith_element, :]
            element_coordinates = self.domain.nodes_coordinates[element_connectivity.astype(int)]
            element_length =  element_coordinates[-1] - element_coordinates[0]
            K_element = np.zeros((self.reference_element.number_nodes,self.reference_element.number_nodes))
            #gauss points loop
            for i in range(self.reference_element.number_gauss_points):
                shape_function_derivative_i = self.reference_element.shape_functions_derivatives_at_gauss_points[i,:] * 2/element_length
                weight_i = self.reference_element.gauss_points_weights[i] * element_length/2
                K_element += weight_i * (shape_function_derivative_i.reshape(-1,1) @ shape_function_derivative_i.reshape(1,-1))
            #assembly
            for e_i, i in zip([0,1], element_connectivity.astype(int)):
                for e_j, j in zip([0,1],element_connectivity.astype(int)):
                    self.K[i,j] +=  K_element[e_i,e_j]

        plt.spy(self.K)
        plt.title(f'Matrix K sparsity pattern', fontsize=15, fontweight='bold')
        plt.show()

    def Solve(self):
        self.d = np.zeros((self.domain.number_of_nodes,self.number_of_steps))
        K = self.K[1:,1:]
        F = np.zeros((self.domain.number_of_nodes -1, 1))
        for i in range(self.number_of_steps):
            F[-1] = i/self.number_of_steps
            self.d[1:,i] = np.squeeze(np.linalg.solve(K,F))

    def SaveResults(self):
        np.save('SnapshotsMatrix.npy', self.d)




class Rom_Simulation(FOM_simulation):

    def __init__(self, basis):
        super().__init__()
        self.basis = basis

    def ComputeSystemMatrix(self):
        self.K = np.zeros((np.shape(self.basis)[1], np.shape(self.basis)[1]))

        #element loop
        for ith_element in range(self.domain.number_of_elements):
            element_connectivity = self.domain.connectivity[ith_element, :]
            element_coordinates = self.domain.nodes_coordinates[element_connectivity.astype(int)]
            element_length =  element_coordinates[-1] - element_coordinates[0]
            K_element = np.zeros((self.reference_element.number_nodes,self.reference_element.number_nodes))
            element_basis = self.basis[element_connectivity.astype(int),:]

            #gauss points loop
            for i in range(self.reference_element.number_gauss_points):
                shape_function_derivative_i = self.reference_element.shape_functions_derivatives_at_gauss_points[i,:] * 2/element_length
                weight_i = self.reference_element.gauss_points_weights[i] * element_length/2
                K_element += weight_i * (shape_function_derivative_i.reshape(-1,1) @ shape_function_derivative_i.reshape(1,-1))

            #assembly
            self.K +=  element_basis.T @ K_element @element_basis

    def Solve(self):
        self.d = np.zeros((self.domain.number_of_nodes,self.number_of_steps))
        F = np.zeros((self.domain.number_of_nodes, 1))

        for i in range(self.number_of_steps):
            F[-1] = i/self.number_of_steps
            F_rom = self.basis.T @ F
            q = np.linalg.solve(self.K,F_rom) #reduced variables
            self.d[:,i] = np.squeeze(self.basis @ q)

    def SaveResults(self):
        np.save('SnapshotsMatrixROM.npy', self.d)

def visualize_resuts_fom(d):
    fig = plt.figure()
    camera = Camera(fig)
    for i in range(np.shape(d)[1]):
        plt.plot(np.linspace(0,np.shape(d)[0],np.shape(d)[0]), d[:,i], 'b-o',alpha=0.5)
        plt.title(f'FOM Simulation Results', fontsize=15, fontweight='bold')
        camera.snap()
    animation = camera.animate()
    animation.save('FOM results.gif', writer = 'imagemagick')
    plt.show()

def visualize_resuts_rom(d, d_rom):
    fig = plt.figure()
    camera = Camera(fig)
    for i in range(np.shape(d_rom)[1]):
        plt.plot(np.linspace(0,np.shape(d)[0],np.shape(d)[0]), d[:,i], 'b-o',alpha=0.5,label='FOM')
        plt.plot(np.linspace(0,np.shape(d_rom)[0],np.shape(d_rom)[0]), d_rom[:,i], 'ro',label='ROM')
        if i==0:
            plt.legend()
        plt.title(f'FOM vs ROM Simulation Results', fontsize=15, fontweight='bold')
        camera.snap()
    animation = camera.animate()
    animation.save('ROM vs FOM results.gif', writer = 'imagemagick')
    plt.show()

def compute_basis(SnapshotsMatrix):
    #taking the svd
    u,s,v = np.linalg.svd(SnapshotsMatrix,full_matrices=False)
    plt.plot(s, 'bo-')
    plt.title(f'Singular Values', fontsize=15, fontweight='bold')
    plt.show()

    u = u[:,0]
    u = u.reshape(-1,1) #since a single mode is obtained, setting expeceted format

    return u


if __name__ == '__main__':
    fom_simulation = FOM_simulation()
    fom_simulation.Run()
    basis = compute_basis(fom_simulation.d)
    rom_simulation = Rom_Simulation(basis)
    rom_simulation.Run()
    visualize_resuts_rom(fom_simulation.d, rom_simulation.d)








