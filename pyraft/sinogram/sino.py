import os
import ctypes
import numpy
import sys
import gc
from ..raftypes import *

import pynfft.nfft as nf
import warnings

##############
# |  pyraft |#
# |  class  |#
##############

# TODO:	Take sinogram's t_extent and image's corners in
#	consideration

def make_fourier_slice_radon_transp( sino, shape = None, sino_zp_factor = 1.5, nfft_zp_factor = 1.2, nfft_m = 2 ):
   """
   make_fourier_slice_radon_transp( shape, sino, sino_zp_factor = 1.5, nfft_zp_factor = 1.2, nfft_m = 2 )

   Creates routines for computing the discrete Radon Transform and its conjugate.

   Parameters:

   sino           : Sinogram. Data is unimportant, but size and ranges are not.
   shape          : Shape of the backprojection result (same of projection argument)
                    You can provide an image instead, in which case its FOV must
                    be the unit square.
   sino_zp_factor : Zero-padding factor for Fourier transform of projection
   nfft_zp_factor : Zero-padding factor for Fourier transform of image
   nfft_m         : Number of summation terms in nfft series.

   Returns:

   return fast_radon, fast_radon_transpose

   functions for computing the Radon transform and its numerical transpose.
  """

   if shape is None:
      shape = ( sino.shape[ 0 ], sino.shape[ 0 ] )
   img = image( shape )

   if ( sino.top_left[ 1 ] != 1.0 ) or \
      ( sino.bottom_right[ 1 ] != -1.0 ):
         raise ValueError( 'Invalid sinogram t range. Must be ( -1.0, 1.0 )' )
   if ( img.top_left != ( -1.0, 1.0 ) ) or \
      ( img.bottom_right != ( 1.0, -1.0 ) ):
         print img.top_left, img.bottom_right
         raise ValueError( 'Invalid image range. Must be ( -1.0, 1.0 ) x ( -1.0, 1.0 ).' )
   if ( img.shape[ 0 ] != sino.shape[ 0 ] ) or \
      ( img.shape[ 1 ] != sino.shape[ 0 ] ):
         warning.warn( 'Attention: img and sino should preferably have matching dimensions.' )

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

   def fast_radon( img ):
      """
         Compute projection through projection-slice theorem
      """

      # Execute plan:
      plan.f_hat = img
      fsino = plan.trafo()

      # Assemble sinogram:
      fsino = numpy.reshape( fsino, ( sino_padded_size, sino.shape[ 1 ] ) )

      # Inverse FFT:
      result = numpy.fft.ifft( fsino, axis = 0 )

      # Shift result:
      result = numpy.fft.ifftshift( result, axes = ( 0, ) )

      # Remove padding:
      result = result[ delta + odd_sino : result.shape[ 0 ] - delta - extra + odd_sino ]

      # Get real part:
      result = numpy.real( result )

      # Normalize:
      result /= ( 0.5 * ( sino.shape[ 0 ] - 1 ) )

      # Return image with appropriate bounding box:
      return image( result, top_left = sino.top_left, bottom_right = sino.bottom_right )

   def fast_radon_transpose( sino ):
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

   return fast_radon, fast_radon_transpose


##############
# |  pyraft |#
# |  class  |#
##############

class radon_method:
   """Methods for projection/backprojection operations."""

   BRESENHAM = 1,
   SLANTSTACK = 2,
   HIERARCHICAL = 3,
   BST = 4,
   ANDERSSON = 5,
   ANDERSSONPB = 6

###############
# |  pyraft  |#
# | function |#
###############


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


###############
# |  pyraft  |#
# | function |#
###############

import pynfft.nfft as nf

def radon_fourier_slice( img, shape=(367, 180), **kwargs ): 

   # Create pyraft.img to hold sinogram:
   sino = image(shape, **kwargs)
   sino.top_left = [0, 1]
   sino.bottom_right = [numpy.pi, -1]
   #SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)

   ####

   sino_zp_factor = 1.5
   nfft_zp_factor = 1.2 
   nfft_m = 2

   # Padded sinogram size (make it even):
   sino_padded_size = int( math.ceil( sino_zp_factor * shape[ 0 ] ) )
   sino_padding_size = sino_padded_size - shape[ 0 ]

   # Fourier radii of FFT irregular samples:
   rhos = numpy.reshape( numpy.fft.fftfreq( sino_padded_size ), ( sino_padded_size, 1 ) )
   # Angular positions of irregular samples:
   thetas = numpy.linspace( sino.top_left[ 0 ], sino.bottom_right[ 0 ], shape[ 1 ] )
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
   Compute projection through projection-slice theorem
   """

   # Execute plan:
   plan.f_hat = img
   fsino = plan.trafo()

   # Assemble sinogram:
   fsino = numpy.reshape( fsino, ( sino_padded_size, sino.shape[ 1 ] ) )

   # Inverse FFT:
   result = numpy.fft.ifft( fsino, axis = 0 )

   # Shift result:
   result = numpy.fft.ifftshift( result, axes = ( 0, ) )

   # Remove padding:
   result = result[ delta + odd_sino : result.shape[ 0 ] - delta - extra + odd_sino ]

   # Get real part:
   result = numpy.real( result )

   # Normalize:
   result /= ( 0.5 * ( sino.shape[ 0 ] - 1 ) )

   # Return image with appropriate bounding box:
   return image( result, top_left = sino.top_left, bottom_right = sino.bottom_right )

	


