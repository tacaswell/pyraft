import os
import ctypes
import numpy 
import sys
from ..raftypes import *
from ..sinogram.sino import *
import gc
import matplotlib.pyplot as plt

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

###############
# |  pyraft  |#
# | function |#
###############

def backprojection_plan( sino_shape=(256,180), shape=(256, 256), padding=2, method=radon_method.BST):
   """Computes the Backprojection using one of several methods"""
  
   # input shape = sinogram shape
   nrays = sino_shape[0]
   nviews= sino_shape[1]

   if method == radon_method.BRESENHAM:
      tupla = []
         
   elif method == radon_method.BST:

      # sinogram @ polar coordinates

      polarsino = image( numpy.zeros([int(nrays/2), 2*nviews]))
      polarsino.top_left=(-numpy.pi, 1.0)
      polarsino.bottom_right=(numpy.pi, 0.0)


      # Zero-padded sinogram

      nrays_zp = libraft.snapTransformSize_bst( ctypes.c_int( padding * int(nrays/2.) - 1) )

      zero_padded_sino = image ( numpy.zeros([nrays_zp, 2*nviews ]) )
    	
      ds = 2.0/float(nrays)
      s_max = nrays_zp * ds
 
      zero_padded_sino.bottom_right = polarsino.bottom_right
      zero_padded_sino.top_left     = ( polarsino.top_left[0], s_max)
 
      # Fourier image of zero padded sinogram (polar representation)

      fftp_re = image( numpy.zeros([nrays_zp, 2*nviews]) )
      fftp_im = image( numpy.zeros([nrays_zp, 2*nviews]) )
  
      fftp_re.bottom_right = zero_padded_sino.bottom_right
      fftp_re.top_left     = zero_padded_sino.top_left

      fftp_im.bottom_right = zero_padded_sino.bottom_right
      fftp_im.top_left     = zero_padded_sino.top_left

      # Fourier image, interpolated to the Cartesian coordinates
  
      fftc_re = image( numpy.zeros( [nrays_zp, nrays_zp] ) )
      fftc_im = image( numpy.zeros( [nrays_zp, nrays_zp] ) )
  
      s1 = max(fftp_re.bottom_right[1], fftp_re.top_left[1])
        
      x0 = -s1
      x1 =  s1
      y0 = -s1
      y1 =  s1

      fftc_re.top_left = [x0, y1]
      fftc_re.bottom_right = [x1, y0]
  
      fftc_im.top_left = [x0, y1]
      fftc_im.bottom_right = [x1, y0]
  
      # Postprocessing and interpolation
      int1 = image( numpy.zeros( [int(fftc_re.shape[0]/2.), int(fftc_re.shape[1]/2.)] ) )
  
      coeff = float(polarsino.shape[0])/float(zero_padded_sino.shape[0])
      
      Cx = coeff * fftc_re.shape[0]
      Cy = coeff * fftc_re.shape[1]
 
      cut  = image( numpy.zeros( [int(Cx), int(Cy)] ) )
        
      int1.top_left = fftc_re.top_left
      int1.bottom_right = fftc_re.bottom_right 

      cut.top_left = ( coeff * fftc_re.top_left[0] , coeff * fftc_re.top_left[1] )
      cut.bottom_right = ( coeff * fftc_re.bottom_right[0], coeff * fftc_re.bottom_right[1] )

      tupla = (method, shape, padding, polarsino, zero_padded_sino, fftp_re, fftp_im, fftc_re, fftc_im, int1, cut, )
 
   elif method == radon_method.ANDERSSON:
     
      # sinogram @ polar coordinates
      nx = shape_res[0]
      ny = shape_res[1]

      prays = int(nrays/2)
      pviews = nviews*2
      polarsino = image( numpy.zeros([prays, pviews]))
      polarsino.top_left=(-numpy.pi, 1.0)
      polarsino.bottom_right=(numpy.pi, 0.0)

      # sinogram @ logpolar coordinates
      ds = 1.0/prays
      dx = 2.0/nx
      dy = 2.0/ny
      dec_st = min(dx,dy)
      r0 = padding*numpy.log(dec_st)

      lrays = libraft.snapTransformSize_bst(int(r0/numpy.log(1-ds)) - 1)
      logpolarsino = image( numpy.zeros( [lrays, pviews] ) ) 
      logpolarsino.top_left = (-numpy.pi, 0.0)
      logpolarsino.bottom_right = (numpy.pi, r0)

      acc = 0.005

      r=numpy.linspace(0, -r0, lrays)
      t=numpy.linspace(-numpy.pi, numpy.pi, pviews)

      [rr,tt] = numpy.meshgrid(r, t)
      z  = numpy.exp(rr)*numpy.cos(tt) - 1
      mask = numpy.abs(z) < acc
      
      kernel = image(numpy.transpose((mask).astype(numpy.double)/acc))
      kernel.top_left = (-numpy.pi,0.0)
      kernel.bottom_right = (numpy.pi,r0)

      # result @ logpolar coordinates
      logpolarres = image( [lrays,pviews] )
      logpolarres.top_left = (-numpy.pi,0)
      logpolarres.bottom_right = (numpy.pi,r0)

      tupla = (method, shape, padding, r0, polarsino, logpolarsino, kernel, logpolarres, )

   elif method == radon_method.ANDERSSONPB:
      tupla = (method, shape, padding, )

   elif method == radon_method.SLANTSTACK:
      tupla = (method, shape, padding, )
   
   else:
      raise TypeError('Unsupported method for Backprojection!')

   return tupla



###############
# |  pyraft  |#
# | function |#
###############

def backprojection(sino, plan, **kwargs):
   """Computes the Backprojection using one of several methods"""

   method = plan[0] 
   shape = plan[1]
   
   # Create pyraft.img to hold backprojection:
   img = image(shape, **kwargs)
   IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)

   # Compute discrete Backprojection:
   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)
   
   if method == radon_method.BRESENHAM:
      libraft.raft_backprojection_bresenham(SINO, IMG, ctypes.c_int(nthreads))
         
   elif method == radon_method.BST:
                
      PLAN = make_RAFT_BST ( plan )

      libraft.raft_backprojection_bst(SINO, IMG, PLAN, ctypes.c_int(nthreads))
      
   elif method == radon_method.ANDERSSON:

      PLAN = make_RAFT_LP( plan ) 

      libraft.raft_backprojection_logpolar(SINO, IMG, PLAN, ctypes.c_int(nthreads) )
      
   elif method == radon_method.ANDERSSONPB:
      libraft.full_sectoral(SINO, IMG , ctypes.c_int(nthreads), ctypes.c_int(4) )	

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


##############
# |  pyraft  |#
# | function |#
##############

def backprojection_gpu(sino, shape=(256, 256), **kwargs):

    img_size = shape[0]

    nrays = sino.shape[0]
    nviews = sino.shape[1]

    sino_buff = numpy.frombuffer(sino.reshape(nrays*nviews)).astype('float32')
    sino_p = sino_buff.ctypes.data_as(POINTER(c_float))

    img = numpy.zeros(img_size*img_size).astype('float32')
    img_p = img.ctypes.data_as(POINTER(c_float))

    libraft.raft_backprojection_slantstack_gpu(img_p, sino_p, img_size, nrays, nviews)

    img_final = img.reshape([img_size,img_size]).astype('float64')

    del img
    gc.collect()
	
    return img_final


##############
# |  pyraft  |#
# | function |#
##############


def back_fourier_slice(sino, shape=(512,512), **kwargs ): 

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


