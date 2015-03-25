import os
import sys
import ctypes
import numpy
from ..raftypes import *

##############
# |  pyraft  |#
# | function |#
##############

def _shepp_logan_desc():
   """Description of Shepp-Logan head phantom."""
   PIover10 = math.pi / 10.0
   desc = numpy.array([ [  1.0, 0.69  , 0.92 , 0.0 , 0.0   , 0.0      ],
                         [ -0.8, 0.6624, 0.874, 0.0 , -0.0184, 0.0      ],
                         [ -0.2, 0.1100, 0.31 , 0.22, 0.0   , -PIover10 ],
                         [ -0.2, 0.1600, 0.41 , -0.22, 0.0   , PIover10 ],
                         [  0.1, 0.2100, 0.25 , 0.0 , 0.35  , 0.0      ],
                         [  0.1, 0.0460, 0.046, 0.0 , 0.1   , 0.0      ],
                         [  0.1, 0.0460, 0.046, 0.0 , -0.1   , 0.0      ],
                         [  0.1, 0.0460, 0.023, -0.08, -0.605 , 0.0      ],
                         [  0.1, 0.0230, 0.023, 0.0 , -0.606 , 0.0      ],
                         [  0.1, 0.0230, 0.046, 0.06, -0.605 , 0.0      ] ]
                     )
   return desc

##############
# |  pyraft  |#
# | function |#
##############

def shepplogan(shape=(256, 256), **kwargs):
   """Creates an image with samples from Golosio's XFCT phantom"""

   # Phantom description:
   desc = _shepp_logan_desc()

   # Extract C-struct from numpy.array:
   DESC = make_RAFT_MATRIX(desc)

   # Create pyraft.image:
   img = image(shape, **kwargs)

   # Extract C-struct from pyraft.image:
   IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)

   # Call libraft function!
   libraft.raft_image_phantom_fromdesc(IMG, DESC);

   # Return image:
   return img

##############
# |  pyraft  |#
# | function |#
##############


def golosio1(shape=(256, 256), **kwargs):
    """Creates an image with samples from Golosio's XFCT phantom"""
    
    # Create pyraft.image:
    img = image(shape, **kwargs)
    
    # Extract C-struct from pyraft.image:
    IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)
    
    # Call libraft function!
    libraft.raft_image_phantom_golosio_I(IMG);
    
    # Return image:
    return img

##############
# |  pyraft  |#
# | function |#
##############
    

def golosio2(shape=(256, 256), **kwargs):
    """Creates an image with samples from Golosio's XFCT phantom"""
    
    # Create pyraft.image:
    img = image(shape, **kwargs)
    
    # Extract C-struct from pyraft.image:
    IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)
    
    # Call libraft function!
    libraft.raft_image_phantom_golosio_II(IMG);
    
    # Return image:
    return img

##############
# |  pyraft  |#
# | function |#
##############


def golosio3(shape=(256, 256), **kwargs):
    """Creates an image with samples from Golosio's XFCT phantom"""
    
    # Create pyraft.image:
    img = image(shape, **kwargs)
    
    # Extract C-struct from pyraft.image:
    IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)
    
    # Call libraft function!
    libraft.raft_image_phantom_golosio_III(IMG);
    
    # Return image:
    return img

##############
# |  pyraft  |#
# | function |#
##############

def gear(shape=(256, 256), **kwargs):
    """Creates an image with samples from Gear phantom"""
    
    # Create pyraft.image:
    img = image(shape, **kwargs)
    
    # Extract C-struct from pyraft.image:
    IMG = make_RAFT_IMAGE(img, img.top_left, img.bottom_right)
    
    # Call libraft function!
    libraft.raft_image_phantom_gear(IMG);
    
    # Return image:
    return img


##############
# |  pyraft  |#
# | function |#
##############

def vol_gear(size=(128)):
 
    """ 3D Gear Phantom """

    g = gear((size, size))
    
    G = numpy.frombuffer(g.data).reshape((size, size))	

    f = numpy.zeros((size, size, size))
    
    _beg = int (15 * size / 32)
    _end = int (17 * size / 32)

    for i in range(_beg, _end):
        f[i, :, :] = g 
    
    return f


