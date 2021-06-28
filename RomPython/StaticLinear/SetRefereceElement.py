import numpy as np

class ReferenceElement():

    def __init__(self, degree):
        self.degree = degree
        self.SetUp()

    def SetUp(self):
        if self.degree == 1:
            self.reference_domain = np.array([-1,1])
            self.number_nodes = 2
            self.number_gauss_points = 2
            self.gauss_points_locations = np.array([-1/np.sqrt(3) , 1/np.sqrt(3)])
            self.gauss_points_weights = np.array([1,1])
            self.shape_functions_at_gauss_points =  np.array([(1-self.gauss_points_locations)/2  , (1+self.gauss_points_locations)/2])
            self.shape_functions_derivatives_at_gauss_points  =  np.array ([[-0.5, 0.5 ],[-0.5, 0.5]])
            self.shape_functions_second_derivatives_at_gauss_points = np.array ([[0, 0 ],[0, 0]])
        else:
            raise NameError(f'reference element of the selected degree {self.degree} not implemented')

def setup_reference_element(degree):
    reference_element_object = ReferenceElement(degree)
    return reference_element_object

if __name__ == '__main__':
    valid_reference_element = ReferenceElement(1)
    print(valid_reference_element)
    invalid_reference_element = ReferenceElement(2) #only linear elements implemented

