import numpy as np

def compute_system_matrix(reference_element, domain):

    K = np.zeros((domain.number_of_nodes, domain.number_of_nodes))

    #element loop
    for ith_element in range(domain.number_of_elements):
        element_connectivity = domain.connectivity[ith_element, :]
        element_coordinates = domain.nodes_coordinates[element_connectivity.astype(int)]
        element_length =  element_coordinates[-1] - element_coordinates[0]
        K_element = np.zeros((reference_element.number_nodes,reference_element.number_nodes))

        #gauss points loop
        for i in range(reference_element.number_gauss_points):
            shape_function_derivative_i = reference_element.shape_functions_derivatives_at_gauss_points[i,:] * 2/element_length
            weight_i = reference_element.gauss_points_weights[i] * element_length/2
            K_element += weight_i * (shape_function_derivative_i.reshape(-1,1) @ shape_function_derivative_i.reshape(1,-1))

        #assembly
        for e_i, i in zip([0,1], element_connectivity.astype(int)):
            for e_j, j in zip([0,1],element_connectivity.astype(int)):
                K[i,j] +=  K_element[e_i,e_j]
    return K