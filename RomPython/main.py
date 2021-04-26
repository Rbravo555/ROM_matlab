import numpy as np
from SetRefereceElement import setup_reference_element
from SetDomain import setup_domain
from ComputeSystemMatrix import compute_system_matrix

def main():
    #set up the reference element
    element_degree = 1
    reference_element = setup_reference_element(element_degree)

    #set up the domain
    number_of_elements = 3
    domain = setup_domain(number_of_elements)

    #build system matrix
    K = compute_system_matrix(reference_element, domain)

    #solve system
    K = K[1:,1:]

    number_of_steps = 10
    u = np.zeros((domain.number_of_nodes,number_of_steps))
    F = np.zeros((domain.number_of_nodes -1, 1))

    for i in range(number_of_steps):
        F[-1] = i/number_of_steps
        u[1:,i] = np.squeeze(np.linalg.solve(K,F))
    print(u)

    #animation of results


if __name__ == '__main__':
    main()