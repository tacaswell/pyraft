from ..raftypes import *
from ..sinogram.sino import *
from ..backprojection.backprojection import *

##############
# |  pyraft  |#
# | function |#
##############

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
    s = numpy.flipud(s)
    
    rays = s.shape[0]
    angles = s.shape[1]
    
    img_size = rays
    
    niter=20
    
    if len(args) > 1:
        img_size = args[1]
    
    
    if len(args) > 2:
        niter = args[2]
    
    sino_buff = numpy.frombuffer(s.reshape(s.shape[0]*s.shape[1])).astype('float32')
    sino_p = sino_buff.ctypes.data_as(POINTER(c_float))
    img = numpy.zeros(img_size*img_size).astype('float32')
    img_p = img.ctypes.data_as(POINTER(c_float))

    libraft.raft_em_gpu(img_p, sino_p, img_size, rays, angles, niter)
    
    img_final = img.reshape([img_size,img_size]).astype('float64')
    return img_final



