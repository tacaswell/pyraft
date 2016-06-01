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

def radon360(phantom, shape=(367, 180)):
   """Computes the Radon transform using old RAFT function"""
  
   # Create pyraft.img to hold backprojection:
   sino = image(shape)
   sino.top_left = [0, 1]
   sino.bottom_right = [numpy.pi, -1.0]
   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)

   PHANTOM = make_RAFT_IMAGE(phantom, phantom.top_left, phantom.bottom_right)

   # Radon transform (slant-stack with pthreads):
   
   libraft.oldraft_radon(SINO, PHANTOM)
   
   return sino


##############
# |  pyraft  |#
# | function |#
##############

def radon(dens, trans, fluor, shape=(367, 180)):
   """Computes the Radon transform using old RAFT function"""
  
   # Create pyraft.img to hold backprojection:

   sino = image(shape)
   sino.top_left = [0, 1]
   sino.bottom_right = [numpy.pi, -1.0]
   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)

   DENS  = make_RAFT_IMAGE(dens, dens.top_left, dens.bottom_right)
   TRANS = make_RAFT_IMAGE(trans, trans.top_left, trans.bottom_right)
   FLUOR = make_RAFT_IMAGE(fluor, fluor.top_left, fluor.bottom_right)

   # XFCT Radon transform (slant-stack with pthreads):
   
   libraft.oldraft_radon_xfct(SINO, DENS, TRANS, FLUOR)
   
   return sino


##############
# |  pyraft  |#
# | function |#
##############

def akt(sino, trans, fluor, shape=(64, 64), niter=4, wtype=1):
   """Computes the Radon transform using old RAFT function"""
  
   # Create pyraft.img to hold backprojection:

   rec = image(shape)
   REC = make_RAFT_IMAGE(rec, rec.top_left, rec.bottom_right)

   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)
   TRANS = make_RAFT_IMAGE(trans, trans.top_left, trans.bottom_right)
   FLUOR = make_RAFT_IMAGE(fluor, fluor.top_left, fluor.bottom_right)

   # AKT iterative method for XFCT inversion:
   
   libraft.oldraft_akt_xfct(REC, SINO, TRANS, FLUOR, ctypes.c_int(niter), ctypes.c_int(wtype) )
   
   return rec


##############
# |  pyraft  |#
# | function |#
##############

def em(sino, trans, fluor, shape=(64, 64), niter=4, wtype=1):
   """Computes the Radon transform using old RAFT function"""
  
   # Create pyraft.img to hold backprojection:

   rec = image(shape)
   REC = make_RAFT_IMAGE(rec, rec.top_left, rec.bottom_right)

   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)
   TRANS = make_RAFT_IMAGE(trans, trans.top_left, trans.bottom_right)
   FLUOR = make_RAFT_IMAGE(fluor, fluor.top_left, fluor.bottom_right)

   # Expect.Maximization for XFCT inversion:
   
   libraft.oldraft_em_xfct(REC, SINO, TRANS, FLUOR, ctypes.c_int(niter), ctypes.c_int(wtype) )
   
   return rec

##############
# |  pyraft  |#
# | function |#
##############

def backprojection360(sino, shape=(256, 256)):
   """Computes the Radon transform using old RAFT function"""
  
   # Create pyraft.img to hold backprojection:
   back = image(shape)
   BACK = make_RAFT_IMAGE(back, back.top_left, back.bottom_right)

   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)

   # Backprojection transform (slant-stack with pthreads):
   
   libraft.oldraft_backprojection(BACK, SINO)
   
   return back

##############
# |  pyraft  |#
# | function |#
##############

def backprojection(sino, trans, fluor, shape=(64, 64)):
   """Computes the Backprojection for XFCT transform """
  
   # Create pyraft.img to hold backprojection:

   rec = image(shape)
   REC = make_RAFT_IMAGE(rec, rec.top_left, rec.bottom_right)

   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)
   TRANS = make_RAFT_IMAGE(trans, trans.top_left, trans.bottom_right)
   FLUOR = make_RAFT_IMAGE(fluor, fluor.top_left, fluor.bottom_right)

   # Expect.Maximization for XFCT inversion:
   
   libraft.oldraft_backprojection_xfct(REC, SINO, TRANS, FLUOR )
   
   return rec


##############
# |  pyraft  |#
# | function |#
##############

def fbp360(sino, shape=(256, 256)):
   """Computes the Radon transform using old RAFT function"""
  
   # Create pyraft.img to hold backprojection:
   back = image(shape)
   BACK = make_RAFT_IMAGE(back, back.top_left, back.bottom_right)

   SINO = make_RAFT_IMAGE(sino, sino.top_left, sino.bottom_right)

   # Backprojection transform (slant-stack with pthreads):
   
   libraft.oldraft_fbp360(BACK, SINO)
   
   return back

