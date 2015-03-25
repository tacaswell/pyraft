import os
import ctypes
import numpy 
import sys
from ..raftypes import *
from ..sinogram.sino import *
import gc

##############
# |  pyraft  |#
# | function |#
##############

'''
def backprojection_with_mesh(sino, N, theta,t):

    R = len(t)
    V = len(theta)
    
    S = numpy.transpose(sino).reshape(V*R,1)       
    S = S.astype(numpy.double)
    B = numpy.ones((N,N)).reshape(N*N,1)
    
    _backprojection( S, B, theta, t, N, R, V, 10)

    return B.reshape((N,N))
'''

##############
# |  pyraft  |#
# | function |#
##############

def backprojection(sino, shape=(256, 256), method=radon_method.BRESENHAM, **kwargs):
   """Computes the Backprojection using one of several methods"""

   # Create pyraft.img to hold backprojection:
   img = image(shape, **kwargs)
   IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)

   # Compute discrete Backprojection:
   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)
   if method == radon_method.BRESENHAM:
      libraft.raft_backprojection_bresenham(SINO, IMG, ctypes.c_int(nthreads))
   elif method == radon_method.SLANTSTACK:
      libraft.raft_backprojection_slantstack(SINO, IMG, ctypes.c_int(nthreads))
   else:
      raise TypeError('Unsupported method for Backprojection!')
   return img

##############
# |  pyraft  |#
# | function |#
##############

def radon_transpose(sino, shape=(256, 256), method=radon_method.BRESENHAM, **kwargs):
   """Computes the transpose of discrete Radon Transform using one of several methods"""

   # Create pyraft.img to hold backprojection:
   img = image(shape, **kwargs)
   IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)

   # Compute transpose:
   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)
   if method == radon_method.BRESENHAM:
      libraft.raft_radon_transpose_bresenham(SINO, IMG, ctypes.c_int(nthreads))
   else:
      raise TypeError('Unsupported method for transpose Radon!')
   return img



def back_gpu(*args):

    if len(args) == 0:
        raise TypeError('Not enough arguments!')        

    s = args[0]
    rays = s.shape[0]
    angles = s.shape[1]
    
    img_size = rays
    
    if len(args) > 1:
        img_size = args[1]
    
    # #

    sino_buff = numpy.frombuffer(s.reshape(s.shape[0]*s.shape[1])).astype('float32')
    sino_p = sino_buff.ctypes.data_as(POINTER(c_float))

    img = numpy.zeros(img_size*img_size).astype('float32')
    img_p = img.ctypes.data_as(POINTER(c_float))

    libraft.raft_backprojection_slantstack_gpu(img_p, sino_p, img_size, rays, angles)

    img_final = img.reshape([img_size,img_size]).astype('float64')

    del img
    gc.collect()
	
    return img_final


