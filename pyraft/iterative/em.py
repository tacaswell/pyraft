from ..raftypes import *
from ..sinogram.sino import *
from ..backprojection.backprojection import *

##############
# |  pyraft |#
# |  class  |#
##############

class em_fourier:

   def __init__( self, sino, shape = None, weight = None, x = None,   **kwargs ):

      # Check if shape is compatible with initial image,
      # if both given:
      if ( not ( shape is None ) ) and \
         ( not ( x is None ) ) and \
         ( shape != x.shape ):
            raise ValueError( 'Incompatible shapes.' )

      # Default shape depends on sinogram spatial resolution:
      if shape is None:
         shape = ( sino.shape[ 0 ], sino.shape[ 0 ] )

      # Initial image:
      if x is None:
         self.x = image( numpy.ones( shape ) )
      else:
         # We already know that x.shape == shape
         self.x = x.copy()

      # If not already an image:
      if not isinstance( self.x, image ):
         self.x = image( self.x )

      # Same for sinogram:
      # TODO: Arredondar o algoritmo e remover listra vertical.
      # TODO: Por enquanto girar por um mínimo irrisório resolve:
      if not isinstance( sino, image ):
         sino = image( sino, x_extent = ( 1e-15, numpy.pi + 1e-15 ) )

      # Sinogram:
      self.sino = sino

      # Create transform functions:
      self.radon, self.radon_transpose = make_fourier_slice_radon_transp( sino, shape, **kwargs )

      if weight is None:
         # Weight will be the application of transpose to vector of ones:
         one = image(
            numpy.ones( sino.shape ),
            top_left = sino.top_left,
            bottom_right = sino.bottom_right
         )

         # Compute weight:
         self.weight = self.radon_transpose( one )
      else:
         self.weight = weight

      # Expose backprojection operator:
      def backprojection( sino ):

	 #from ..misc import imagesc
	 #mask = ( numpy.abs(self.weight) < 1e-4 )
	 #imagesc(mask)

         tmp = self.radon_transpose( sino )
         tmp /= self.weight

         return tmp

      self.backprojection = backprojection

   def __call__( self, x = None, return_intermediate = False ):

      # No current iteration passed, using internal:

      if x is None:
         x = self.x

      # Makes sure we work with an image:
      if not isinstance( x, image ):
         x = image( x.copy() )

      # Compute Radon transform:
      tp_i = time.time()
      Ax = self.radon( x )
      tp_f = time.time() - tp_i
      
      # Iteration:
      tr_i = time.time()
      bp = self.backprojection( self.sino / Ax )
      x *= bp
      tr_f = time.time() - tp_i
      
      if return_intermediate:
         # TODO: Is Ax.copy() required, or will Ax do?
         return x.copy(), Ax.copy()#, tp_f, tr_f
      else:
         return x.copy(), Ax.copy()#, tp_f, tr_f


###############
# |  pyraft  |#
# | function |#
###############

def make_em(sino, shape, projection_method=radon_method.BRESENHAM, backprojection_method=radon_method.BRESENHAM, weight=None, **kwargs):
   
   # We have to compute weight:
   if not weight:
      one = image(
         numpy.ones(sino.shape),
         top_left=sino.top_left,
         bottom_right=sino.bottom_right
      )

      # If shape argument is an image, we need to create
      # a weight with same geometry:
      if isinstance(shape, image):
         weight = radon_transpose(
            one,
            shape.shape,
            top_left=shape.top_left,
            bottom_right=shape.bottom_right,
            method=backprojection_method,
            **kwargs
         )
      # Otherwise we use the provided arguments:
      else:
         weight = radon_transpose(one, shape, method=backprojection_method, **kwargs)

   # EM algorithm:
   def em(x):
      # Compute Radon transform:
      Ax = radon(x, sino.shape, top_left=sino.top_left, bottom_right=sino.bottom_right, method=projection_method)
      # Iteration:
      Ax = sino / Ax
      x = x * radon_transpose(Ax, x.shape, top_left=x.top_left, bottom_right=x.bottom_right, method=backprojection_method) / weight

      return x

   return em

###############
# |  pyraft  |#
# | function |#
###############

def em(sino, shape=None, niter=5, **kwargs):

   if shape == None:
      shape = [sino.shape[1], sino.shape[1]]

   if isinstance(shape, image):
   # Initial image was provided:
      if kwargs:
         raise TypeError('Unhandled arguments!')
      x = shape
   else:
   # Use random starting image:
      x = image(shape, **kwargs)
      x[:] = 1.0

   _em_iter = make_em(sino, x)
   for i in range(0, niter): 
      x = _em_iter(x)
   return x

###############
# |  pyraft  |#
# | function |#
###############

def em_gpu(*args):

    if len(args) == 0:
        raise TypeError('Not enough arguments!')        

    s = args[0]
    
    rays = s.shape[0]
    angles = s.shape[1]
    
    img_size = rays
    
    niter = 20
    
    if len(args) > 1:
        img_size = args[1]
    
    
    if len(args) > 2:
        niter = args[2]
    
    sino_buff = numpy.frombuffer(s.reshape(s.shape[0]*s.shape[1])).astype('float32')
    sino_p = sino_buff.ctypes.data_as(POINTER(c_float))
    img = numpy.zeros(img_size*img_size).astype('float32')
    img_p = img.ctypes.data_as(POINTER(c_float))

    libraft.raft_em_gpu(img_p, sino_p, img_size, rays, angles, niter)
    
    img_final = img.reshape([img_size,img_size])

    return img_final


###############
# |  pyraft  |#
# | function |#
###############

#
# OS-EM using Fourier-Slice for Radon and Backprojection transforms
# Author(s): 	Camila de Lima
# 		Elias S. Helou
#

def osem(sino, shape=None, niter=5, nblocks=4, **kwargs):
  
   '''
   if shape == None:
      shape = [sino.shape[0], sino.shape[0]] 

   if isinstance(shape, image):
      # Initial image was provided:
      if kwargs:
         raise TypeError('Unhandled arguments!')
      x = shape
   else:
      # Use random starting image:
      x = image(shape, **kwargs)
      x[:] = 1.0
   '''
	
   x = image ( numpy.ones( (sino.shape[0], sino.shape[0])) )
   
   b = 0
   e = (sino.shape[1])/nblocks
   e0 = e
   dim_i = 0
   dim_f = (numpy.math.pi)/nblocks
   dim_f0 = dim_f 

   sinograma = [None]*nblocks
   radon = [None]*nblocks
   radontransp = [None]*nblocks
   alg = [None]*nblocks

   for i in range(nblocks):
      
      sinograma[i] = image( sino[:,b:e] , x_extent = ( 1e-15 + dim_i, 1e-15 + dim_f))
      #sinograma[i].top_left = ( dim_i, 1)
      #sinograma[i].bottom_right = (-dim_f, -1)      

      radon[i], radontransp[i] = make_fourier_slice_radon_transp(sinograma[i])
   
      alg[i] = em_fourier(sinograma[i])

      #from ..misc import imagesc
      #imagesc(radontransp[i](sinograma[i]))

      b += e0
      e += e0
      dim_i += dim_f0
      dim_f += dim_f0
      
   #

   for k in range (niter):
      for i in range (nblocks):
         x, Ax = alg[i](x)
	
   return x


