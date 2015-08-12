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
         
   elif method == radon_method.BST:
      libraft.raft_backprojection_miqueles(SINO, IMG , ctypes.c_int(nthreads))
   
   elif method == radon_method.BST:
      libraft.raft_backprojection_andersson(SINO, IMG , ctypes.c_int(nthreads))
   
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


##############
# |  pyraft  |#
# | function |#
##############


def back_fourier_slice( sino, shape=(512,512), **kwargs ): 

   # Create pyraft.img to hold backprojection:
   img = image(shape, **kwargs)
 
   ####

   sino_zp_factor = 1.5
   nfft_zp_factor = 1.2 
   nfft_m = 2

   # Padded sinogram size (make it even):
   sino_padded_size = int( math.ceil( sino_zp_factor * sino.shape[ 0 ] ) )
   sino_padding_size = sino_padded_size - sino.shape[ 0 ]

   # Fourier radii of FFT irregular samples:
   rhos = numpy.reshape( numpy.fft.fftfreq( sino_padded_size ), ( sino_padded_size, 1 ) )
   # Angular positions of irregular samples:
   thetas = numpy.linspace( sino.top_left[ 0 ], sino.bottom_right[ 0 ], sino.shape[ 1 ] )
   # Trigonometric values:
   trig_vals = numpy.reshape(
      numpy.transpose(
         numpy.array( [ numpy.sin( thetas ), -numpy.cos( thetas ) ] )
         ),
      ( 1, thetas.shape[ 0 ] * 2 )
      )
   # Finally, desired irregular samples:
   sample_positions = numpy.reshape( rhos * trig_vals, ( thetas.shape[ 0 ] * rhos.shape[ 0 ], 2 ) )

   # Computations later required to remove padding.
   delta = sino_padding_size / 2
   extra = sino_padding_size % 2
   odd_sino = sino.shape[ 0 ] % 2

   # Plan nonuniform FFT:
   plan = nf.NFFT(
      d = 2,
      N = img.shape,
      M = sample_positions.shape[ 0 ],
      flags = ( 'PRE_PHI_HUT', 'PRE_PSI' ),
      n = [ i * nfft_zp_factor for i in img.shape  ],
      m = nfft_m
      )
   plan.x = sample_positions
   plan.precompute()

   """
   Compute backprojection through projection-slice theorem
   """

   # Zero-pad sinogram:
   fsino = numpy.zeros( ( sino_padded_size, sino.shape[ 1 ] ) )
   fsino[ delta + odd_sino : fsino.shape[ 0 ] - delta - extra + odd_sino ] = sino

   # Shift sinogram columns:
   fsino = numpy.fft.fftshift( fsino, axes = ( 0, ) )

   # Fourier transform of projections:
   fsino = numpy.fft.fft( fsino, axis = 0 )

   # Dissasemble transformed sinogram:
   plan.f = numpy.reshape( fsino, ( fsino.shape[ 0 ] * fsino.shape[ 1 ], 1 ) )

   # Compute adjoint of nouniform Fourier transform:
   result = plan.adjoint()

   # Get real part:
   result = numpy.real( result )

   # Normalize:
   result /= ( 0.5 * ( sino_padded_size - 1 ) * ( img.shape[ 1 ] - 1 ) )

   return image( result, top_left = img.top_left, bottom_right = img.bottom_right )


