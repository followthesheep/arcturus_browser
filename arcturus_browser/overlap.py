import numpy as np

def overlap(a, b):
    # return the indices in a that overlap with b, also returns the corresponding index in b
    # only works if both a and b are unique! This is not very efficient but it works
    bool_a = np.in1d(a,b)
    ind_a = np.arange(len(a))
    ind_a = ind_a[bool_a]

    ind_b = np.array([np.argwhere(b == a[x]) for x in ind_a]).flatten()
    return ind_a,ind_b
        

    
