import numpy
from scipy.ndimage.interpolation import rotate, shift
from ..raftypes import image

def center_sino(sino):
    #input sinogram: (rays x angles)	

    N = sino.shape[1]
    R = sino.shape[0]

    th = numpy.linspace(0, numpy.pi, N)
    th.shape = [len(th), 1]
    t = numpy.linspace(-1, 1, R)
    t.shape = [len(t), 1]
    t0 = float(t[0])	
    dt = float(t[1] - t0) 

    a1 = numpy.ones([N, 1])
    a2 = numpy.cos(th)
    a3 = numpy.sin(th)

    A = numpy.hstack([a1, a2, a3])

    T = numpy.kron(t, numpy.ones([1, N]))

    s1 = T * sino
    s = s1.sum(0) / sino.sum(0)
    s.shape = [len(s), 1]

    c = float( numpy.linalg.lstsq(A, s)[0][0] )

    ind_zero = int ( numpy.ceil( (0-t0)/dt ) )
    ind_c = int ( numpy.ceil( (c-t0)/dt ) )

    offset = int(ind_zero - ind_c)

    return numpy.roll( sino, offset, 0) , offset   

