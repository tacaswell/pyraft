import os
import ctypes
import numpy
import sys
import gc
from ..raftypes import *

#############
# |  pyraft |#
# |  class  |#
#############

class radon_method:
   """Methods for projection/backprojection operations."""

   BRESENHAM = 1,
   SLANTSTACK = 2,
   HIERARCHICAL = 3

##############
# |  pyraft  |#
# | function |#
##############


def radon(img, shape=(367, 180), method=radon_method.BRESENHAM, **kwargs):
   """Computes the Radon transform using one of several methods"""

   # Test for Radon-space extents:
   if (not ('top_left' in kwargs)) and \
      (not ('bottom_right' in kwargs)) and \
      (not ('x_extent' in kwargs)) and \
      (not ('y_extent' in kwargs)):
      if 'theta_extent' in kwargs:
         kwargs[ 'x_extent' ] = kwargs[ 'theta_extent' ]
         del kwargs[ 'theta_extent' ]
      else:
         kwargs[ 'x_extent' ] = (0.0, math.pi)
      if 't_extent' in kwargs:
         kwargs[ 'y_extent' ] = kwargs[ 't_extent' ]
         del kwargs[ 't_extent' ]
   # Create pyraft.img to hold sinogram:
   sino = image(shape, **kwargs)
   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)

   # img is a pyraft.image: compute discrete Radon transform:
   if isinstance(img, image):
      IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)
      if method == radon_method.BRESENHAM:
         libraft.raft_radon_bresenham(IMG, SINO, ctypes.c_int(nthreads))
      elif method == radon_method.SLANTSTACK:
         libraft.raft_radon_slantstack(IMG, SINO, ctypes.c_int(nthreads))
      else:
         raise TypeError('Unsupported method for Radon Transform!')
      return sino

   # img is a numpy.array: consider it a description and compute
   # exact Radon transform:
   elif isinstance(img, numpy.ndarray):
      DESC = make_RAFT_MATRIX(img)
      libraft.raft_radon_fromdesc(SINO, DESC)
      return sino


###############
# |  pyraft  |#
# | function |#
###############

def radon_view(*args):

    if len(args) == 0:
        raise TypeError('Not enough arguments!')
            
    f = args[0]
    theta = args[1]

    tmin = f.top_left[0]
    tmax = f.bottom_right[0]
	
    t = numpy.linspace(tmin, tmax, f.shape[1])
    
    if len(args) > 2:
        t = args[2]
                
    view = numpy.zeros(t.shape)

    VIEW = make_RAFT_VECTOR(view)	
    T = make_RAFT_VECTOR(t)	
    F = make_RAFT_IMAGE(f)	
    
    libraft.raft_radon_slantstack_view(VIEW, F, T, theta)
	   
    return view

###############
# |  pyraft  |#
# | function |#
###############

def radon_withmesh(*args):

    if len(args) == 0:
        raise TypeError('Not enough arguments!')
            
    f = args[0]
    theta = args[1]

    tmin = f.top_left[0]
    tmax = f.bottom_right[0]

    t = numpy.linspace(tmin, tmax, f.shape[1])
    
    if len(args) > 2:
        t = args[2]

    sino = image((len(t), len(theta)))
    SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)
                
    THETA = make_RAFT_VECTOR(theta)	
    T = make_RAFT_VECTOR(t)	
    F = make_RAFT_IMAGE(f)	
    
    libraft.raft_radon_slantstack_withmesh(F, SINO, T, THETA, ctypes.c_int(nthreads))
	   
    return sino


###############
# |  pyraft  |#
# | function |#
###############

def radon_3d_stack(*args):

    if len(args) == 0:
        raise TypeError('Not enough arguments!')        

    f = args[0]

    theta = numpy.linspace(0, numpy.pi, 180)
    t = numpy.linspace(-1, 1, f.shape[1])
    
    if len(args) > 1:
        theta = args[1]
    
    if len(args) > 2:
        t = args[2]
    
    # #    

    V = len(theta)
    R = len(t)
    N = f.shape[0]
    S = numpy.zeros((V, N, R))

    for i in range(N):
        slice_ = image(f[:, i, :])
        r = radon_withmesh(slice_ , theta, t)
        S[:, i, :] = numpy.transpose(r)
    return S

###############
# |  pyraft  |#
# | function |#
###############

def snapshot(*args):
    
    if len(args) == 0:
        raise TypeError('Not enough arguments!')            

    f = numpy.asarray(args[0])
    theta = 0
    t = numpy.linspace(-1.0, 1.0, f.shape[1])

    if len(args) == 2:
        theta = float(args[1])    
    
    if len(args) > 2:
        t = args[2]
        
    R = len(t)
    N = f.shape[0]

    snap = numpy.zeros((N, N))
    
    for k in range(N):
        sview = radon_view(image(f[:, k, :]) , theta, t)
        snap[k, :] = sview
        
    return snap
        

###############
# |  pyraft  |#
# | function |#
###############

def radon_gpu(*args):

    if len(args) == 0:
        raise TypeError('Not enough arguments!')        

    f = args[0]
    rays = f.shape[0]
    angles = 360

    if len(args) > 1:
        rays = args[1]
    
    if len(args) > 2:
        angles = args[2]
    
    img_size = f.shape[0]

    f_buff = numpy.frombuffer(f.data).astype('float32')
    f_p = f_buff.ctypes.data_as(POINTER(c_float))
    
    
    ############################################################
    ##radon
    	
    sino = numpy.zeros(rays * angles).astype('float32')
    sino_p = sino.ctypes.data_as(POINTER(c_float))
        
    libraft.raft_radon_slantstack_gpu(sino_p, f_p, img_size, rays, angles)
        
    sino_final2 = sino.reshape([rays, angles])
    sino_final = image(sino_final2.astype('float64'))
    sino_final.top_left = (0,1.0)
    sino_final.bottom_right = (numpy.pi,-1.0)

    del sino    
    gc.collect()
    
    return sino_final


###############
# |  pyraft  |#
# | function |#
###############

def sinogram_cutborder(sino, percentage):
    #input sinogram: (rays x angles)    

    R = int(sino.shape[0])
    
    t = numpy.linspace(-1, 1, R)
    t0 = float(t[0])	
    tN = float(t[R-1])
    dt = float(t[1] - t0) 	
	
    beta = float(percentage/100.00)    

    k = int ( numpy.ceil( (tN-beta-t0)/dt) )

    return sino[R-k:k,:]

