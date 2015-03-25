import os
import numpy 
from .ringfun import *
from ..raftypes import *

##############
# |  pyraft  |#
# | function |#
##############

''' Filtered sinogram for FBP '''

def lowpassino_fbp(*args):
       
    if len(args) == 0:
        print ("error!")
        return 0

    sino = args[0]
            
    tmin = sino.bottom_right[1]
    tmax = sino.top_left[1]

    t = numpy.linspace(tmin, tmax, sino.shape[0])
    
    if len(args) > 1:
        t = args[1]

    R = int(sino.shape[0])
    V = int(sino.shape[1])
    
    dt = (t[R - 1] - t[0]) / (R - 1)
    wc = 1 / (2 * dt)
    w = numpy.linspace(-wc, wc, R)
    h = numpy.abs(w) * numpy.cos(numpy.pi * w * dt)
    
    G = numpy.fft.fftshift(numpy.transpose(numpy.kron(numpy.ones((V, 1)), h)))
    B = numpy.fft.fft(sino, axis=0) 
    
    C = B * G

    #---	
    
    D = numpy.fft.ifft(C, axis=0) 
    
    W = image(D.real)	
    W.top_left = sino.top_left
    W.bottom_right = sino.bottom_right

    return W


##############
# |  pyraft  |#
# | function |#
##############

''' Ring Filters without blocks'''

def ring(SINO):

	d1 = ring_nb(SINO, 1, 1)
	
	d2 = ring_nb(SINO, 2, 1)
	
	p = d1 * d2

	alpha = 1.5

	d = numpy.sqrt(p + alpha * numpy.abs(p.min()))
		
	return d

##############
# |  pyraft  |#
# | function |#
##############

''' Ring Filters with blocks'''

def ring_block(*args):

	SINO = args[0]

	if len(args) == 1:
		print ("No blocks were given, using 2")
		nblocks = 2
	else:
		nblocks = args[1]
   
	size = int(SINO.shape[0] / nblocks)
	
	d1 = ring_b(SINO, 1, 1, size)

	d2 = ring_b(SINO, 2, 1, size)

	p = d1 * d2

	alpha = 1.5
	
	d = numpy.sqrt(p + alpha * numpy.fabs(p.min()))
	
	return d
