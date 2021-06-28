import numpy as np
from matplotlib import pyplot as plt
from celluloid import Camera




def set_up_Newmark_coefficients(dt, alpha=0.25, beta=0.5):
    a0=1/(alpha*(dt**2))
    a1=beta/(alpha*dt)
    a2=1/(alpha*dt)
    a3=(1/(2*alpha))-1
    a4=(beta/alpha)-1
    a5=(dt/2)*((beta/alpha)-2)
    a6=dt*(1-beta)
    a7=beta*dt

    return a0,a1,a2,a3,a4,a5,a6,a7




def Nonlinear_1D_Bar_FOM(number_of_elements,total_time,time_step_size):
    number_of_time_steps =int(total_time/time_step_size)
    a0,a1,a2,a3,a4,a5,a6,a7= set_up_Newmark_coefficients(time_step_size)
    restrained_dofs=[0,] #clamped bc
    deltaX= np.zeros((number_of_elements+1,1))
    tol=1e-6 #Setting tolerance for the nonlinear iterations
    U = np.zeros((number_of_elements+1,1)) #initializing displacements
    U[-1] = 0.1
    Un = U.copy()
    Ud=np.zeros((number_of_elements+1,1))#initializing velocities
    Udn=Ud.copy()
    Udd=np.zeros((number_of_elements+1,1)) #initializing acceleration
    Uddn=Udd.copy()
    F= np.zeros((len(U), 1)) #creating the force of zeros
    #F[number_of_elements]=0.1 #Applying a force at the end of the bar
    K_static_el=(np.array([[2,0],[0,2]]) * float(number_of_elements))  #For nonlinear case f=u**2 -1
    M_el= np.array([[2,1],[1,2]]) / (6. * number_of_elements )
    SnapshotsMatrixDisplacements = np.zeros(((len(U)), number_of_time_steps))
    SnapshotsMatrixVelocities = np.zeros(((len(U)), number_of_time_steps))
    SnapshotsMatrixAccelerations = np.zeros(((len(U)), number_of_time_steps))

    #Time step loop
    for time_step in range (number_of_time_steps):
        print("time step", time_step)

        U=Un.copy()
        Ud=Udn.copy()
        Udd=Uddn.copy()
        Un=U+(Ud*time_step_size)

        Uddn=a0*(Un-U)-a2*Ud-a3*Udd #obtain acceleration using Newmark coefficients
        Udn=Ud+a6*Udd+a7*Uddn #obtain velocity using Newmark coefficients

        deltaX[number_of_elements-1]=50 #Setting a high value of dx to enter loop
        iter=0

        while np.linalg.norm(deltaX)>tol:
            iter+=1
            print("iteration", iter)

            #Creating the global K and R for Newton-Raphson
            Residual_global = np.zeros((number_of_elements+1,1))
            K_dynamic_global = np.zeros((number_of_elements+1,number_of_elements+1))

            #Assemble contributions for every element
            for i in range(number_of_elements):
                F_ext_el=F[i:i+2,0]
                Un_el=Un[i:i+2,0]
                Uddn_el=Uddn[i:i+2,0]

                F_int_el = np.power(Un_el,2)-1 #This is the nonlinear term f(x)

                Residual_element=F_ext_el-F_int_el-(M_el@Uddn_el)

                #Tangent matrix
                K_dynam_el=K_static_el+(a0*M_el)

                auxiliary_variable_assemeble_K=np.zeros((number_of_elements+1,number_of_elements+1))
                auxiliary_variable_assemble_R=np.zeros((number_of_elements+1,1))

                auxiliary_variable_assemeble_K[i:i+2,i:i+2] = K_dynam_el
                auxiliary_variable_assemble_R[i:i+2,0] = Residual_element

                K_dynamic_global+=auxiliary_variable_assemeble_K
                Residual_global+=auxiliary_variable_assemble_R


            #Remove fixed degrees of freedom
            for dof in restrained_dofs:
                K_dynamic_global = np.delete(K_dynamic_global, dof, axis=0)
                K_dynamic_global = np.delete(K_dynamic_global, dof, axis=1)
                Residual_global= np.delete(Residual_global, dof, axis=0)

            #Solve global system
            dx=np.linalg.solve(K_dynamic_global, Residual_global)

            #Correcting displacement
            deltaX[1:number_of_elements+1,0] = dx.transpose()
            Un = Un + deltaX


            Uddn=a0*(Un-U)-a2*Ud-a3*Udd #update acceleration using Newmark coefficients
            Udn=Ud+a6*Udd+a7*Uddn #update velocity using Newmark coefficients

        SnapshotsMatrixDisplacements[:,time_step]=Un.T
        SnapshotsMatrixVelocities[:,time_step]=Udn.T
        SnapshotsMatrixAccelerations[:,time_step]=Uddn.T


    return SnapshotsMatrixDisplacements, SnapshotsMatrixVelocities, SnapshotsMatrixAccelerations






def Nonlinear_1D_Bar_ROM(number_of_elements, total_time, time_step_size, Phi):
    number_of_time_steps =int(total_time/time_step_size)
    number_of_dofs = int (Phi.shape[1])
    a0,a1,a2,a3,a4,a5,a6,a7= set_up_Newmark_coefficients(time_step_size)
    restrained_dofs=[0,] #clamped bc
    tol=1e-6 #Setting tolerance for the nonlinear iterations
    deltaQ= np.zeros((number_of_dofs , 1))
    U = np.zeros((number_of_elements+1,1)) #initializing displacements
    U[-1] = 0.1
    q= Phi.T @ U #initializing displacement
    qn=q.copy()
    qd=np.zeros((number_of_dofs , 1))#initializing velocity
    qdn=qd.copy()
    qdd=np.zeros((number_of_dofs , 1)) #initializing acceleration
    qddn=qdd.copy()
    F= np.zeros((Phi.shape[0], 1)) #creating the force of zeros
    #F[number_of_elements]=0.1 #Applying a force at the end of the bar
    K_static_el=(np.array([[2,0],[0,2]]) * float(number_of_elements))  #For nonlinear case f=u**2 -1
    M_el= np.array([[2,1],[1,2]]) / (6. * number_of_elements )
    SnapshotsMatrixDisplacementsReduced = np.zeros((Phi.shape[1], number_of_time_steps))
    SnapshotsMatrixVelocitiesReduced = np.zeros((Phi.shape[1], number_of_time_steps))
    SnapshotsMatrixAccelerationsReduced = np.zeros((Phi.shape[1], number_of_time_steps))

    #Time step loop
    for time_step in range (number_of_time_steps):
        print("time step", time_step)

        q=qn
        qd=qdn
        qdd=qddn
        qn=q+(qd*time_step_size)

        qddn=a0*(qn-q)-a2*qd-a3*qdd #obtain acceleration using Newmark coefficients
        qdn=qd+a6*qdd+a7*qddn #obtain velocity using Newmark coefficients

        deltaQ[(number_of_dofs-1),0]=50 #Setting a high value of dx to enter loop

        iter=0
        while abs(sum(deltaQ))>tol:
            iter+=1
            print("iteration", iter)

            #Creating the global K and R for NR
            Residual_global = np.zeros((number_of_dofs,1))
            K_dynamic_global = np.zeros((number_of_dofs,number_of_dofs))

            #Assemble contributions for every element
            for i in range(number_of_elements):
                Phi_elem = Phi[i:i+2,:]
                F_ext_el = (F[i:i+2,0:1])
                Un_el = Phi_elem @ qn
                Uddn_el = Phi_elem @ qddn

                #Calculate Residual
                F_int_el=np.power(Un_el,2)-1 #This is the nonlinear term f(x)

                #F_int_el=np.exp(Un_el)
                Residual_element = F_ext_el-F_int_el-(M_el @ Uddn_el)
                Residual_element_svd = Phi_elem.transpose() @ Residual_element


                #Calculate Tangent
                K_dynam_el=K_static_el+(a0*M_el)
                K_dynam_el_rom = Phi_elem.T @ K_dynam_el @ Phi_elem

                K_dynamic_global+=K_dynam_el_rom
                Residual_global+=Residual_element_svd


            #Remove fixed degrees of freedom
            #NO NEED TO REMOVE THEM WHEN USING HOMOGENEOUS BCs

            #Solve global system
            dx = np.linalg.solve(K_dynamic_global, Residual_global)

            #Correct displacement
            deltaQ=dx
            qn=qn + deltaQ

            qddn=a0*(qn-q)-a2*qd-a3*qdd #update acceleration using Newmark coefficients
            qdn=qd+a6*qdd+a7*qddn #update velocity using Newmark coefficients

        SnapshotsMatrixDisplacementsReduced[:,time_step] = qn.T
        SnapshotsMatrixVelocitiesReduced[:,time_step] = qdn.T
        SnapshotsMatrixAccelerationsReduced[:,time_step] = qddn.T

    SnapshotsMatrixDisplacements = Phi @ SnapshotsMatrixDisplacementsReduced
    SnapshotsMatrixVelocities = Phi @ SnapshotsMatrixVelocitiesReduced
    SnapshotsMatrixAccelerations = Phi @ SnapshotsMatrixAccelerationsReduced

    return SnapshotsMatrixDisplacements, SnapshotsMatrixVelocities, SnapshotsMatrixAccelerations






def plotting_solution(dis_FOM, dis_ROM, vel_FOM, vel_ROM, acc_FOM, acc_ROM, node_id):

    plt.figure()
    plt.plot(dis_FOM[node_id,:], 'go', label='FOM')
    plt.plot(dis_ROM[node_id,:], 'r', label='ROM')
    plt.title('Displacement comparison')

    plt.figure()
    plt.plot(vel_FOM[node_id,:], 'go', label='FOM')
    plt.plot(vel_ROM[node_id,:], 'r', label='ROM')
    plt.title('Velocity comparison')

    plt.figure()
    plt.plot(acc_FOM[node_id,:], 'go', label='FOM')
    plt.plot(acc_ROM[node_id,:], 'r', label='ROM')
    plt.title('Acceleration comparison')

    plt.show()






def compute_percentual_difference(Original, Approximation):
    return (np.linalg.norm(Original - Approximation) / np.linalg.norm(Original)) * 100




def compute_basis(SnapshotsMatrix, truncation_tolerance):
    #taking the svd
    u,s,v = np.linalg.svd(SnapshotsMatrix,full_matrices=False)

    #truncating matrix of left singular vectors
    DOWN =np.sum(np.power(s,2))
    UP=DOWN.copy()
    for i in range(len(s)):
        UP = UP - s[i]**2
        if np.sqrt(UP/DOWN)<truncation_tolerance:
            Phi = u[:,:i]
            break

    return Phi




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


if __name__=='__main__':

    # TODO there are some interesting stuff happening
    # Is the artificial nonlinear term making the example unstable
    # depending on the number of elements, the stability lasts more or less
    # 20 elements and 1.5 secs give a stable solution

    #simulation_parameters
    number_of_elements = 10
    total_time = 1.2
    time_step_size = 0.01

    #run FOM
    dis_FOM,vel_FOM,acc_FOM = Nonlinear_1D_Bar_FOM(number_of_elements,total_time,time_step_size)

    # #see animation of FOM solution
    # visualize_resuts_fom(dis_FOM)

    #build basis
    truncation_tolerance = 5e-2
    Phi = compute_basis(dis_FOM, truncation_tolerance)

    #run ROM
    dis_ROM, vel_ROM, acc_ROM = Nonlinear_1D_Bar_ROM(number_of_elements,total_time,time_step_size, Phi)

    #plot solution at a node
    node_id = -1 #the node at the free end of the bar
    plotting_solution(dis_FOM, dis_ROM, vel_FOM, vel_ROM, acc_FOM, acc_ROM, node_id)

    #check approximation
    print( 'displacement difference: ', compute_percentual_difference(dis_FOM, dis_ROM), '%')
    print( 'velocity difference: ', compute_percentual_difference(vel_FOM,vel_ROM) , '%')
    print( 'acceleration difference: ', compute_percentual_difference(acc_FOM,acc_ROM), '%')


    #see animation of FOM solution
    visualize_resuts_rom(dis_FOM, dis_ROM)


