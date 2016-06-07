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


###############
# |  pyraft  |#
# | function |#
###############

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

	


