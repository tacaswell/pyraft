#We use the one encoding: utf8
import ctypes
from ctypes import *
import ctypes.util
import multiprocessing
import math
import os
import sys
import numpy
import time

nthreads = multiprocessing.cpu_count()

# Load required libraies:

#XFCT_LIBS


libstdcpp = ctypes.CDLL( ctypes.util.find_library( "stdc++" ), mode=ctypes.RTLD_GLOBAL )
libblas   = ctypes.CDLL( ctypes.util.find_library( "blas" ), mode=ctypes.RTLD_GLOBAL )
libfftw3  = ctypes.CDLL( ctypes.util.find_library( "fftw3" ), mode=ctypes.RTLD_GLOBAL )
libfftw3_threads  = ctypes.CDLL( ctypes.util.find_library( "fftw3_threads" ), mode=ctypes.RTLD_GLOBAL )

#############
#|  pyraft |#
#|  import |#
#############

_lib = "lib/libraft"

if sys.version_info[0] >= 3:
    import sysconfig
    ext = sysconfig.get_config_var('SO')
else:
    ext = '.so'

_path = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + _lib + ext
libraft  = ctypes.CDLL(_path)


#############
#|  pyraft |#
#|  class  |#
#############

# Define a special type for the 'double *' argument
class DoubleArrayType:
    def from_param(self, param):
        typename = type(param).__name__
        if hasattr(self, 'from_' + typename):
            return getattr(self, 'from_' + typename)(param)
        elif isinstance(param, ctypes.Array):
            return param
        else:
            raise TypeError("Can't convert %s" % typename)
            
    # Cast from array.array objects
    def from_array(self, param):
        if param.typecode != 'd':
            raise TypeError('must be an array of doubles')
            ptr, _ = param.buffer_info()
            return ctypes.cast(ptr, ctypes.POINTER(ctypes.c_double))
            
    # Cast from lists/tuples
    def from_list(self, param):
        val = ((ctypes.c_double)*len(param))(*param)
        return val

    from_tuple = from_list
    
    def from_ndarray(self, param): # Cast from a numpy array
        return param.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


DoubleArray = DoubleArrayType()

#############
#|  pyraft |#
#|  class  |#
#############

# "double *" type:
_c_double_p = ctypes.POINTER(ctypes.c_double)

class RAFT_MATRIX( ctypes.Structure ):
   """A raft_matrix from raft:"""
   _fields_ = [ ( "p_data", _c_double_p ),
                ( "lines", ctypes.c_int ),
                ( "line_stride", ctypes.c_int ),
                ( "columns", ctypes.c_int ),
                ( "column_stride", ctypes.c_int )
               ]

class RAFT_IMAGE( ctypes.Structure ):
   """A raft_image from raft:"""
   _fields_ = [ ( "data", RAFT_MATRIX ),
                ( "tl_x", ctypes.c_double ),
                ( "tl_y", ctypes.c_double ),
                ( "br_x", ctypes.c_double ),( "br_y", ctypes.c_double )
               ]

class RAFT_VECTOR ( ctypes.Structure ):
    """A raft_vector from raft:"""
    _fields_ = [ ( "p_data", _c_double_p ),
                 ( "size", ctypes.c_int ),
                 ( "stride", ctypes.c_int )               
               ]
    
class RAFT_BST ( ctypes.Structure ):
    """A raft_bst plan from raft:"""
    _fields_ = [( "padding_coeff", ctypes.c_double ),
		( "polarsino", RAFT_IMAGE ),
		( "zero_padded_sino", RAFT_IMAGE ),
		( "fftp_re", RAFT_IMAGE ),
		( "fftp_im", RAFT_IMAGE ),
		( "fftc_re", RAFT_IMAGE ),
		( "fftc_im", RAFT_IMAGE ),
		( "int1", RAFT_IMAGE ),
		( "cut", RAFT_IMAGE )
               ]


class RAFT_LP ( ctypes.Structure ):
    """A raft_plan_logpolaar plan from raft:"""
    _fields_ = [( "padding_coeff", ctypes.c_double ),
                ( "r0", ctypes.c_double ),
		( "polar_sino", RAFT_IMAGE ),
		( "logpolar_sino", RAFT_IMAGE ),
		( "kernel", RAFT_IMAGE ),
		( "logpolar_res", RAFT_IMAGE )
               ]

#########################
#|        pyraft       |#
#| Function prototypes |#
#########################

libraft.raft_image_phantom_fromdesc.argtypes = [ RAFT_IMAGE, RAFT_MATRIX ]
libraft.raft_image_phantom_fromdesc.restype = None
libraft.raft_radon_fromdesc.argtypes = [ RAFT_IMAGE, RAFT_MATRIX ]
libraft.raft_radon_fromdesc.restype = None

libraft.raft_image_phantom_gear.argtypes = [ RAFT_IMAGE ]
libraft.raft_image_phantom_gear.restype = None
libraft.raft_image_phantom_golosio_I.argtypes = [ RAFT_IMAGE ]
libraft.raft_image_phantom_golosio_I.restype = None
libraft.raft_image_phantom_golosio_II.argtypes = [ RAFT_IMAGE ]
libraft.raft_image_phantom_golosio_II.restype = None
libraft.raft_image_phantom_golosio_III.argtypes = [ RAFT_IMAGE ]
libraft.raft_image_phantom_golosio_III.restype = None

libraft.raft_radon_bresenham.argtypes = [ RAFT_IMAGE, RAFT_IMAGE, ctypes.c_int ]
libraft.raft_radon_bresenham.restype = None
libraft.raft_backprojection_bresenham.argtypes = [ RAFT_IMAGE, RAFT_IMAGE, ctypes.c_int ]
libraft.raft_backprojection_bresenham.restype = None
libraft.raft_radon_transpose_bresenham.argtypes = [ RAFT_IMAGE, RAFT_IMAGE, ctypes.c_int ]
libraft.raft_radon_transpose_bresenham.restype = None

libraft.raft_radon_slantstack.argtypes = [ RAFT_IMAGE, RAFT_IMAGE, ctypes.c_int ]
libraft.raft_radon_slantstack.restype = None
libraft.raft_backprojection_slantstack.argtypes = [ RAFT_IMAGE, RAFT_IMAGE, ctypes.c_int ]
libraft.raft_backprojection_slantstack.restype = None

libraft.raft_backprojection_bst.argtypes = [ RAFT_IMAGE, RAFT_IMAGE, RAFT_BST, ctypes.c_int ]
libraft.raft_backprojection_bst.restype = None

libraft.raft_backprojection_logpolar.argtypes = [ RAFT_IMAGE, RAFT_IMAGE, RAFT_LP, ctypes.c_int ]
libraft.raft_backprojection_logpolar.restype = None

libraft.raft_haar.argtypes = [ RAFT_MATRIX, ctypes.c_int, ctypes.c_int ]
libraft.raft_haar.restype = None

libraft.raft_radon_slantstack_view.argtypes = [RAFT_VECTOR, RAFT_IMAGE, RAFT_VECTOR, ctypes.c_double]
libraft.raft_radon_slantstack_view.restype = None

libraft.raft_radon_slantstack_withmesh.argtypes = [RAFT_IMAGE, RAFT_IMAGE, RAFT_VECTOR, RAFT_VECTOR, ctypes.c_int]
libraft.raft_radon_slantstack_withmesh.restype = None

libraft.raft_radon_slantstack_withmesh.argtypes = [RAFT_IMAGE, RAFT_IMAGE, RAFT_VECTOR, RAFT_VECTOR, ctypes.c_int]
libraft.raft_radon_slantstack_withmesh.restype = None

libraft.raft_radon_slantstack_withmesh.argtypes = [RAFT_IMAGE, RAFT_IMAGE, RAFT_VECTOR, RAFT_VECTOR, ctypes.c_int]
libraft.raft_radon_slantstack_withmesh.restype = None

libraft.raft_filter1d_ramp.argtypes = [ ctypes.c_double, RAFT_IMAGE, ctypes.c_int ]
libraft.raft_filter1d_ramp.restype = None

libraft.snapTransformSize_bst.argtypes = [ ctypes.c_int ]
libraft.snapTransformSize_bst.restype = ctypes.c_int 

libraft.raft_kernel_lp_create.argtypes = [ RAFT_IMAGE, ctypes.c_double ]
libraft.raft_kernel_lp_create.restype = None

libraft.raft_image_s2p.argtypes = [RAFT_IMAGE, RAFT_IMAGE, ctypes.c_int ]
libraft.raft_image_s2p.restypes = None

libraft.raft_image_p2c.argtypes = [RAFT_IMAGE, RAFT_IMAGE, ctypes.c_int ]
libraft.raft_image_p2c.restypes = None

#
#CUDA PART
#

try:
    libraft.raft_backprojection_slantstack_gpu.argtypes
    libraft.raft_backprojection_slantstack_gpu.argtypes = [POINTER(c_float), POINTER(c_float), c_size_t, c_size_t, c_size_t]   
    libraft.raft_backprojection_slantstack_gpu.restype= None   
    libraft.raft_radon_slantstack_gpu.argtypes = [POINTER(c_float), POINTER(c_float), c_size_t, c_size_t, c_size_t]   
    libraft.raft_em_gpu.restype= None  
    libraft.raft_radon_slantstack_gpu.restype= None  
    libraft.raft_em_gpu.argtypes = [POINTER(c_float), POINTER(c_float), c_size_t, c_size_t, c_size_t, c_size_t]   
    libraft.raft_tr_em_gpu.argtypes = [POINTER(c_float), POINTER(c_float), POINTER(c_float), c_size_t, c_size_t, c_size_t, c_size_t]   
    libraft.raft_tr_em_gpu.restype= None 
except:
    print ('PyRAfT: No CUDA Functions!!')


##########################################################################
#
#			XFCT FUNCTIONS ( with GSL & CONFUSE )
#
##########################################################################

try:
    libraft.oldraft_radon.argtypes = [RAFT_IMAGE, RAFT_IMAGE ]
    libraft.oldraft_radon.restypes = None
    
    libraft.oldraft_radon_xfct.argtypes = [RAFT_IMAGE, RAFT_IMAGE, RAFT_IMAGE, RAFT_IMAGE ]
    libraft.oldraft_radon_xfct.restypes = None
    
    libraft.oldraft_akt_xfct.argtypes = [RAFT_IMAGE, RAFT_IMAGE, RAFT_IMAGE, RAFT_IMAGE, ctypes.c_int, ctypes.c_int ]
    libraft.oldraft_akt_xfct.restypes = None
    
    libraft.oldraft_em_xfct.argtypes = [RAFT_IMAGE, RAFT_IMAGE, RAFT_IMAGE, RAFT_IMAGE, ctypes.c_int, ctypes.c_int ]
    libraft.oldraft_em_xfct.restypes = None
    
    libraft.oldraft_backprojection.argtypes = [RAFT_IMAGE, RAFT_IMAGE]
    libraft.oldraft_backprojection.restypes = None
    
    libraft.oldraft_backprojection_xfct.argtypes = [RAFT_IMAGE, RAFT_IMAGE, RAFT_IMAGE, RAFT_IMAGE]
    libraft.oldraft_backprojection_xfct.restypes = None
    
    libraft.oldraft_fbp360.argtypes = [RAFT_IMAGE, RAFT_IMAGE]
    libraft.oldraft_fbp360.restypes = None
except:
    print ('PyRAfT: No XFCT Functions!!')
    

#############
#|  pyraft |#
#|  class  |#
#############

class image( numpy.ndarray ):
   """This class represents an image. This may not be the ideal foundation, 
   because there are already some options for image classes. Further study is necessary."""

   def __new__(
                subtype,
                shape,
                top_left = None,
                bottom_right = None,
                extent = None,
                x_extent = None,
                y_extent = None,
                **kwargs
              ):
      """Creates and returns a new object of the correct subtype"""

      # Which field-of-view arguments where given?
      extent_given = extent or x_extent or y_extent
      corner_given = top_left or bottom_right

      # Are they acceptable?
      if extent_given and corner_given:
         raise TypeError( 'Mutually exclusive arguments given.' )

      # Extent given, adjust corners:
      if extent_given:
         # Extent given by parts:
         if not extent:
            if not x_extent:
               x_extent = ( -1.0, 1.0 )
            if not y_extent:
               y_extent = ( -1.0, 1.0 )
         # Extent given fully:
         else:
            x_extent = ( extent[ 0 ], extent[ 1 ] )
            y_extent = ( extent[ 2 ], extent[ 3 ] )
         # Finally, we can set up corners:
         top_left     = ( x_extent[ 0 ], y_extent[ 1 ] )
         bottom_right = ( x_extent[ 1 ], y_extent[ 0 ] )

      # pyraft.image given as argument
      if isinstance( shape, image ):

         # Check for given corners:
         if not extent_given:
            if not top_left:
               top_left = shape.top_left
            if not bottom_right:
               bottom_right = shape.bottom_right

         # No arguments other than corners can be taken:
         if kwargs:
            raise TypeError( 'Unhandled arguments!' )

         ## In here, shape is actually a pyraft.image:
         #obj = numpy.asarray( shape ).view( subtype )
         # TODO: No view, make a copy! But there must be a neater way...
         obj = numpy.ndarray.__new__( subtype, shape.shape, **kwargs )
         obj[ ... ] = shape[ ... ]

      else:

         # Check for given corners:
         if not extent_given:
            if not top_left:
               top_left = ( -1.0, 1.0 )
            if not bottom_right:
               bottom_right = ( 1.0, -1.0 )

         # numpy.ndarray given as argument:
         if isinstance( shape, numpy.ndarray ):

            if kwargs:
            # No arguments other than corners can be taken:
               raise TypeError( 'Unhandled arguments!' )

            # In here, shape is actually a numpy.ndarray:
            #obj = numpy.asarray( shape ).view( subtype )
            # TODO: No view, make a copy! But there must be a neater way...
            obj = numpy.ndarray.__new__( subtype, shape.shape, **kwargs )
            obj[ ... ] = shape[ ... ]

         # We must create a zero array:
         else:

            # Default data type is double:
            if not ( 'dtype' in kwargs ):
               kwargs[ 'dtype' ] = numpy.float64
            obj = numpy.ndarray.__new__( subtype, shape, **kwargs )
            obj[ : ] = 0.0

      # All relevant dimensions must match:
      if ( len( obj.shape ) != len( top_left ) ) or ( len( top_left ) != len( bottom_right ) ):
         raise TypeError( 'Dimensions must match!' )

      # Set new attributes:
      obj.top_left = top_left
      obj.bottom_right = bottom_right
      try:
         obj.sampling_distances = ( ( bottom_right[ 0 ] - top_left[ 0 ] ) / ( obj.shape[ 1 ] - 1.0 ),
                                    ( bottom_right[ 1 ] - top_left[ 1 ] ) / ( obj.shape[ 0 ] - 1.0 )
                                  )
      except ZeroDivisionError:
         obj.sampling_distances = ( 0.0, 0.0 )
      return obj

   def __array_finalize__( self, obj ):
      """Set self attributes"""
      if obj is None: return # When ran from __new__

      # Else do the job:
      self.top_left = getattr( obj, 'top_left', None )
      self.bottom_right = getattr( obj, 'bottom_right', None )
      self.sampling_distances = getattr( obj, 'sampling_distances', None )

   def __reduce__( self ):

      # Initial state is only ndarray state:
      full_state = list( numpy.ndarray.__reduce__( self ) )

      #Further attributes:
      image_state = ( self.top_left, self.bottom_right, self.sampling_distances )

      # Add image attributes:
      full_state[ 2 ] = ( full_state[ 2 ], image_state )

      return tuple( full_state )

   def __setstate__( self, state ):

      # Call superclass' __setstate__:
      numpy.ndarray.__setstate__( self, state[ 0 ] )

      # Set our own state:
      self.top_left = state[ 1 ][ 0 ]
      self.bottom_right = state[ 1 ][ 1 ]
      self.sampling_distances = state[ 1 ][ 2 ]

   def sample_coordinates( self, idx ):
      """ Returns coordinates of sample """
      return ( self.top_left[ 0 ] + idx[ 1 ] * self.sampling_distances[ 0 ], self.top_left[ 1 ] + idx[ 0 ] * self.sampling_distances[ 1 ] )

   def get_y_coordinate( self, idx ):
      """ Returns y-coordinate of row """
      return self.top_left[ 1 ] + idx * self.sampling_distances[ 1 ]

   def get_x_coordinate( self, idx ):
      """ Returns x-coordinate of column """
      return self.top_left[ 0 ] + idx * self.sampling_distances[ 0 ]

   # Extent:
   def extent( self ):
      return ( self.top_left[ 0 ], self.bottom_right[ 0 ], self.bottom_right[ 1 ], self.top_left[ 1 ] )

##############
#|  pyraft  |#
#| function |#
##############

def make_RAFT_MATRIX( array ):
   """Make a raft_matrix from a numpy.ndarray"""
   return RAFT_MATRIX( ctypes.cast( array.ctypes.data, _c_double_p ), 
                       ctypes.c_int( int(array.shape[ 0 ] ) ),
                       ctypes.c_int( int(array.strides[ 1 ] / 8 ) ),
                       ctypes.c_int( int(array.shape[ 1 ] ) ),
                       ctypes.c_int( int(array.strides[ 0 ] / 8 ) )
                      )

##############
#|  pyraft  |#
#| function |#
##############

def make_RAFT_IMAGE( array, top_left = ( -1.0, 1.0 ), bottom_right = ( 1.0, -1.0 ) ):
   """Make a raft_matrix from a numpy.ndarray from a pyraft.image or from a pyraft.RAFT_MATRIX"""
   if isinstance( array, numpy.ndarray ):
      return RAFT_IMAGE( make_RAFT_MATRIX( array ),
                         ctypes.c_double( top_left[ 0 ] ),
                         ctypes.c_double( top_left[ 1 ] ),
                         ctypes.c_double( bottom_right[ 0 ] ),
                         ctypes.c_double( bottom_right[ 1 ] )
                       )
   elif isinstance( array, RAFT_MATRIX ):
      return RAFT_IMAGE( array,
                         ctypes.c_double( top_left[ 0 ] ),
                         ctypes.c_double( top_left[ 1 ] ),
                         ctypes.c_double( bottom_right[ 0 ] ),
                         ctypes.c_double( bottom_right[ 1 ] )
                       )
   elif isinstance( array, image ):
      return RAFT_IMAGE( array,
                         ctypes.c_double( array.top_left[ 0 ] ),
                         ctypes.c_double( array.top_left[ 1 ] ),
                         ctypes.c_double( array.bottom_right[ 0 ] ),
                         ctypes.c_double( array.bottom_right[ 1 ] )
                       )

##############
#|  pyraft  |#
#| function |#
##############

def make_RAFT_BST( tupla ):
   """Make a raft_bst from a tuple of {numpy.ndarray, pyraft.image or pyraft.matrix} """
   
   return RAFT_BST(   ctypes.c_double(tupla[2]), 
		      make_RAFT_IMAGE(tupla[3]), 
		      make_RAFT_IMAGE(tupla[4]), 
		      make_RAFT_IMAGE(tupla[5]), 
		      make_RAFT_IMAGE(tupla[6]), 
		      make_RAFT_IMAGE(tupla[7]), 
		      make_RAFT_IMAGE(tupla[8]), 
		      make_RAFT_IMAGE(tupla[9]), 
		      make_RAFT_IMAGE(tupla[10]) 
		    )

def make_RAFT_LP( tupla ):
   """Make a raft_plan_logpolar from a tuple of {numpy.ndarray, pyraft.image or pyraft.matrix} """
   
   return RAFT_LP(  ctypes.c_double(tupla[2]),
                    ctypes.c_double(tupla[3]),
                    make_RAFT_IMAGE(tupla[4]),
                    make_RAFT_IMAGE(tupla[5]),
                    make_RAFT_IMAGE(tupla[6]),
                    make_RAFT_IMAGE(tupla[7])
           )



##############
#|  pyraft  |#
#| function |#
##############

def make_RAFT_VECTOR( array ):
   """Make a raft_vector from a numpy.ndarray"""
   return RAFT_VECTOR( ctypes.cast( array.ctypes.data, _c_double_p ), 
                       ctypes.c_int( int(array.shape[ 0 ]) ),
                       ctypes.c_int( int(array.strides[ 0 ] / 8) )
                     )
   
##############
#|  pyraft  |#
#| function |#
##############

def img2array( raftimage ):
    
    total = (raftimage.data.lines * raftimage.data.columns)

    local = ctypes.cast(raftimage.data.p_data, ctypes.POINTER(ctypes.c_double*total))

    return numpy.frombuffer(local.contents).reshape( (raftimage.data.lines, raftimage.data.columns) )
	

def vec2array( raftvector ):

    local = ctypes.cast(raftvector.p_data, ctypes.POINTER(ctypes.c_double*raftvector.size))

    return numpy.frombuffer(local.contents)


##############
#|          |#
##############

if __name__ == "__main__":
   pass

