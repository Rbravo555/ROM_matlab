import numpy as np

class Domain():
    #1D domain class
    def __init__(self, number_of_elements):
        self.number_of_elements = number_of_elements
        self.SetUp()


    def SetUp(self):
        self.number_of_nodes = self.number_of_elements+1
        self.nodes_coordinates = np.linspace(0,1,self.number_of_nodes)
        self.connectivity = np.c_[np.linspace(0,self.number_of_nodes-2,self.number_of_nodes-1), np.linspace(1,self.number_of_nodes-1,self.number_of_nodes-1) ]


def setup_domain(number_of_elements):
    domain_object = Domain(number_of_elements)
    return domain_object





