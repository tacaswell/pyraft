import numpy
from ..backprojection.backprojection import *
from ..filters.filters import *
from ..raftypes import img2array

##############
# |  pyraft  |#
# | function |#
##############

def fbp(*args):
    
    if len(args) == 0:
        print ("pyraft:: filtered backprojection Error")
        return 1
        
    sino = args[0]
    N = sino.shape[0]
    
    filt_sino = lowpassino_fbp(sino)

    img = backprojection(filt_sino , shape=(N, N))

    return img 


##############
# |  pyraft  |#
# | function |#
##############

def fbp_gpu(*args):

    if len(args) == 0:
        raise TypeError('Not enough arguments!')        

    s = args[0]
    rays = s.shape[0]
    angles = s.shape[1]
    
    img_size = rays
    
    if len(args) > 1:
        img_size = args[1]
    
    ##

    filtsino = lowpassino_fbp(s)
    filtsino_buff = numpy.frombuffer(filtsino.reshape(filtsino.shape[0]*filtsino.shape[1])).astype('float32')
    filtsino_p =  filtsino_buff.ctypes.data_as(POINTER(c_float))

    img = numpy.zeros(img_size*img_size).astype('float32')
    img_p = img.ctypes.data_as(POINTER(c_float))

    t = time.time()
    libraft.raft_backprojection_slantstack_gpu(img_p, filtsino_p, img_size, rays, angles)
    
    #img_final = img.reshape([img_size,img_size]).astype('float64')
    
    return img_final

'''
def fbp_withmesh(sino, N, theta,t):

    fsino = lowpassino_fbp(sino)
    
    return backprojection_with_mesh (fsino, theta, t)
'''

