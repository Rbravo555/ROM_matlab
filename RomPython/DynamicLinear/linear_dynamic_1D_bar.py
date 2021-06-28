import numpy as np
from matplotlib import pyplot as plt
from celluloid import Camera

from scipy.linalg import eigh





class Domain(object):
    #1D domain class
    def __init__(self, number_of_elements):
        self.number_of_elements = number_of_elements
        self.SetUp()

    def SetUp(self):
        self.number_of_nodes = self.number_of_elements+1
        self.nodes_coordinates = np.linspace(0,1,self.number_of_nodes)
        self.connectivity = np.c_[np.linspace(0,self.number_of_nodes-2,self.number_of_nodes-1), np.linspace(1,self.number_of_nodes-1,self.number_of_nodes-1) ]





class ReferenceElement(object):

    def __init__(self):
        self.SetUp()

    def SetUp(self):
        self.reference_domain = np.array([-1,1])
        self.number_nodes = 2
        self.number_gauss_points = 2
        self.gauss_points_locations = np.array([-1/np.sqrt(3) , 1/np.sqrt(3)])
        self.gauss_points_weights = np.array([1,1])
        self.shape_functions_at_gauss_points =  np.array([(1-self.gauss_points_locations)/2  , (1+self.gauss_points_locations)/2])
        self.shape_functions_derivatives_at_gauss_points  =  np.array ([[-0.5, 0.5 ],[-0.5, 0.5]])
        self.shape_functions_second_derivatives_at_gauss_points = np.array ([[0, 0 ],[0, 0]])




class FOM_simulation(object):

    def __init__( self, number_of_elements=5, total_time=20,time_step_size=.1):
        self.total_time = total_time
        self.number_of_time_steps =int(total_time/time_step_size)
        self.time_step_size = time_step_size
        self.setup_domain(number_of_elements)
        self.setup_reference_element()


    def setup_reference_element(self):
        self.reference_element = ReferenceElement()


    def setup_domain(self, number_of_elements):
        self.domain = Domain(number_of_elements)


    def Run(self):
        self.ComputeSystemMatrix()
        self.Solve()


    def ComputeSystemMatrix(self):
        self.K = np.zeros((self.domain.number_of_nodes, self.domain.number_of_nodes))
        self.M = self.K.copy()
        #element loop
        for ith_element in range(self.domain.number_of_elements):
            element_connectivity = self.domain.connectivity[ith_element, :]
            element_coordinates = self.domain.nodes_coordinates[element_connectivity.astype(int)]
            element_length =  element_coordinates[-1] - element_coordinates[0]
            K_element = np.zeros((self.reference_element.number_nodes,self.reference_element.number_nodes))
            M_element = K_element.copy()
            #gauss points loop
            for i in range(self.reference_element.number_gauss_points):
                shape_function_i = self.reference_element.shape_functions_at_gauss_points[i,:]
                shape_function_derivative_i = self.reference_element.shape_functions_derivatives_at_gauss_points[i,:] * 2/element_length
                weight_i = self.reference_element.gauss_points_weights[i] * element_length/2
                K_element += weight_i * (shape_function_derivative_i.reshape(-1,1) @ shape_function_derivative_i.reshape(1,-1))
                M_element += weight_i * (shape_function_i.reshape(-1,1) @ shape_function_i.reshape(1,-1))
            #assembly
            for e_i, i in zip([0,1], element_connectivity.astype(int)):
                for e_j, j in zip([0,1],element_connectivity.astype(int)):
                    self.M[i,j] +=  M_element[e_i,e_j]
                    self.K[i,j] +=  K_element[e_i,e_j]

        plt.spy(self.K)
        plt.title(f'Matrix K sparsity pattern', fontsize=15, fontweight='bold')
        plt.figure()
        plt.spy(self.M)
        plt.title(f'Matrix M sparsity pattern', fontsize=15, fontweight='bold')
        plt.show()


    def GetInitialDisplacement(self, applied_force):
        applied_foce_vector = np.zeros((self.domain.number_of_elements))
        applied_foce_vector[-1] = applied_force
        Displacement_old = np.zeros((self.domain.number_of_nodes))
        Displacement_old[1:] = np.squeeze(np.linalg.solve(self.K[1:,1:], applied_foce_vector))

        return Displacement_old


    def Solve(self):
        self.set_up_Newmark_coefficients()

        applied_force = 1
        Displacement_old = self.GetInitialDisplacement(applied_force) #solve static problem to determine original deformation state

        Displacement_new = Displacement_old.copy()
        Velocity_old = np.zeros((self.domain.number_of_nodes))
        Acceleration_old = np.zeros((self.domain.number_of_nodes))

        self.SnapshotsMatrixDisplacements = np.zeros((np.shape(Displacement_old)[0],self.number_of_time_steps))
        self.SnapshotsMatrixDisplacements[:,0] = Displacement_old
        self.SnapshotsMatrixVelocities = np.zeros((np.shape(Displacement_old)[0],self.number_of_time_steps))
        self.SnapshotsMatrixVelocities[:,0] = Velocity_old
        self.SnapshotsMatrixAccelerations = np.zeros((np.shape(Displacement_old)[0],self.number_of_time_steps))
        self.SnapshotsMatrixAccelerations[:,0] = Acceleration_old

        F = np.zeros((self.domain.number_of_nodes))
        K_hat = self.K + self.a0*self.M

        for i in range(1,self.number_of_time_steps):
            F_hat = F + self.M @ (self.a0*Displacement_old  +  self.a2*Velocity_old   + self.a3 * Acceleration_old)

            #solve system
            Displacement_new[1:] = np.linalg.solve(K_hat[1:,1:], F_hat[1:])

            #update acceleration using Newmark coefficients
            Acceleration_new = self.a0 * (Displacement_new - Displacement_old)  - self.a2*Velocity_old -  self.a3 * Acceleration_old

            #update velocity using Newmark coefficients
            Velocity_new = Velocity_old + self.a6 * Acceleration_old + self.a7*Acceleration_new

            self.SnapshotsMatrixDisplacements[:,i] = Displacement_new
            self.SnapshotsMatrixVelocities[:,i] = Velocity_new
            self.SnapshotsMatrixAccelerations[:,i] = Acceleration_new

            #reset_variables
            Displacement_old = Displacement_new.copy()
            Velocity_old = Velocity_new.copy()
            Acceleration_old = Acceleration_new.copy()


    def set_up_Newmark_coefficients(self, alpha=0.25, beta=0.5):
        dt = self.time_step_size
        self.a0=1/(alpha*(dt**2))
        self.a1=beta/(alpha*dt)
        self.a2=1/(alpha*dt)
        self.a3=(1/(2*alpha))-1
        self.a4=(beta/alpha)-1
        self.a5=(dt/2)*((beta/alpha)-2)
        self.a6=dt*(1-beta)
        self.a7=beta*dt



class ROM_simulation(FOM_simulation):


    def __init__(self, number_of_elements,total_time,time_step_size,basis):
        super().__init__(number_of_elements, total_time, time_step_size)
        self.basis = basis

    def ComputeSystemMatrix(self):
        self.K = np.zeros((np.shape(self.basis)[1], np.shape(self.basis)[1]))
        self.M = self.K.copy()
        #element loop
        for ith_element in range(self.domain.number_of_elements):
            element_connectivity = self.domain.connectivity[ith_element, :]
            element_coordinates = self.domain.nodes_coordinates[element_connectivity.astype(int)]
            element_length =  element_coordinates[-1] - element_coordinates[0]
            K_element = np.zeros((self.reference_element.number_nodes,self.reference_element.number_nodes))
            M_element = K_element.copy()
            element_basis = self.basis[element_connectivity.astype(int),:]

            #gauss points loop
            for i in range(self.reference_element.number_gauss_points):
                shape_function_i = self.reference_element.shape_functions_at_gauss_points[i,:]
                shape_function_derivative_i = self.reference_element.shape_functions_derivatives_at_gauss_points[i,:] * 2/element_length
                weight_i = self.reference_element.gauss_points_weights[i] * element_length/2
                K_element += weight_i * (shape_function_derivative_i.reshape(-1,1) @ shape_function_derivative_i.reshape(1,-1))
                M_element += weight_i * (shape_function_i.reshape(-1,1) @ shape_function_i.reshape(1,-1))
            #assembly
            self.K +=  element_basis.T @ K_element @element_basis
            self.M +=  element_basis.T @ M_element @element_basis


    def GetInitialDisplacement(self, applied_force):
        applied_foce_vector = np.zeros((self.domain.number_of_nodes))
        applied_foce_vector[-1] = applied_force
        applied_foce_vector_rom = self.basis.T @ applied_foce_vector

        Displacement_old_rom = np.squeeze(np.linalg.solve(self.K , applied_foce_vector_rom))

        return Displacement_old_rom


    def Solve(self):
        self.set_up_Newmark_coefficients()

        applied_force = 1
        Displacement_old_rom = self.GetInitialDisplacement(applied_force) #solve static problem to determine original deformation state

        Displacement_new_rom = Displacement_old_rom.copy()
        Velocity_old_rom = np.zeros((np.shape(self.basis)[1]))
        Acceleration_old_rom = np.zeros((np.shape(self.basis)[1]))

        self.SnapshotsMatrixDisplacements_rom = np.zeros((np.shape(self.basis)[1],self.number_of_time_steps))
        self.SnapshotsMatrixDisplacements_rom[:,0] = Displacement_old_rom
        self.SnapshotsMatrixVelocities_rom = np.zeros((np.shape(self.basis)[1],self.number_of_time_steps))
        self.SnapshotsMatrixVelocities_rom[:,0] = Velocity_old_rom
        self.SnapshotsMatrixAccelerations_rom = np.zeros((np.shape(self.basis)[1],self.number_of_time_steps))
        self.SnapshotsMatrixAccelerations_rom[:,0] = Acceleration_old_rom

        F = np.zeros((np.shape(self.basis)[1]))
        K_hat = self.K + self.a0*self.M

        for i in range(1,self.number_of_time_steps):
            F_hat = F + self.M @ (self.a0*Displacement_old_rom  +  self.a2*Velocity_old_rom  + self.a3 * Acceleration_old_rom)

            #solve system
            Displacement_new_rom  = np.linalg.solve(K_hat, F_hat)

            #update acceleration using Newmark coefficients
            Acceleration_new_rom = self.a0 * (Displacement_new_rom - Displacement_old_rom)  - self.a2*Velocity_old_rom -  self.a3 * Acceleration_old_rom

            #update velocity using Newmark coefficients
            Velocity_new_rom = Velocity_old_rom + self.a6 * Acceleration_old_rom + self.a7*Acceleration_new_rom

            self.SnapshotsMatrixDisplacements_rom[:,i] = Displacement_new_rom
            self.SnapshotsMatrixVelocities_rom[:,i] = Velocity_new_rom
            self.SnapshotsMatrixAccelerations_rom[:,i] = Acceleration_new_rom

            #reset_variables
            Displacement_old_rom = Displacement_new_rom.copy()
            Velocity_old_rom = Velocity_new_rom.copy()
            Acceleration_old_rom = Acceleration_new_rom.copy()






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




def compute_basis(SnapshotsMatrix, truncation_tolerance=1e-4):
    #taking the svd
    u,s,v = np.linalg.svd(SnapshotsMatrix,full_matrices=False)
    norm_S = np.linalg.norm(SnapshotsMatrix)

    Phi=None

    #truncating matrix of left singular vectors
    UP = np.sum(np.power(s,2))
    DOWN= 0

    for i in range(len(s)):
        DOWN += s[i]**2
        UP -= s[i]**2
        if np.sqrt(UP/DOWN) < truncation_tolerance*norm_S:
            Phi = u[:,:i+1]
            break

    if Phi is None:
        Phi = u

    return Phi


if __name__ == '__main__':
    #simulations parameters
    number_of_elements = 4
    total_time = 4
    time_step_size = .1


    #run FOM simulation
    fom_simulation = FOM_simulation(number_of_elements,total_time,time_step_size)
    fom_simulation.Run()


    #POD basis
    #truncation_tolerance_for_svd = 1e-4
    #basis = compute_basis(fom_simulation.SnapshotsMatrixDisplacements, truncation_tolerance_for_svd)

    #Modal analysis basis
    _, basis =  eigh(fom_simulation.K[1:,1:],fom_simulation.M[1:,1:])
    basis = np.r_[np.zeros((1,np.shape(basis)[1])), basis]

    #run ROM simulation
    rom_simulation = ROM_simulation(number_of_elements,total_time,time_step_size,basis)
    rom_simulation.Run()

    ##visualize the displacements  (SLOW )
    visualize_resuts_rom(fom_simulation.SnapshotsMatrixDisplacements, rom_simulation.basis @  rom_simulation.SnapshotsMatrixDisplacements_rom )

    print('The difference is FOM vs ROM is : ', np.linalg.norm(fom_simulation.SnapshotsMatrixDisplacements - rom_simulation.basis @  rom_simulation.SnapshotsMatrixDisplacements_rom) / np.linalg.norm(fom_simulation.SnapshotsMatrixDisplacements), '%')








