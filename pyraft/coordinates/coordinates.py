import os
import sys
import ctypes
import numpy
from ..raftypes import *

##############
# |  pyraft  |#
# | function |#
##############

def sinogram_to_polar(sino):
   """Transforming a sinogram to polar coordinates"""

   R = sino.shape[0]
   V = sino.shape[1]

   # Sinogram @ polar coordinates:
   img = image([R/2, 2*V], top_left=[0,1], bottom_right=[2*numpy.pi,0])
   IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)

   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)
   

   # Call libraft function!
   libraft.raft_image_s2p(SINO, IMG, ctypes.c_int(nthreads) );

   # Return image:
   return img



def polar_to_cartesian(N, img):
   """Transform a polar image from polar to cartesian coordinates"""

   # image @ cartesian coordinates:
   cart = image( [N, N] )
   CART = make_RAFT_IMAGE(cart, cart.top_left, cart.bottom_right)

   IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)
   
   # Call libraft function!
   libraft.raft_image_p2c(IMG, CART, ctypes.c_int(nthreads) );

   # Return image:
   return cart
