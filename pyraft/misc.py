import os
import numpy 
from pylab import *
import matplotlib.pyplot as plt

################################################ 
#               Python codes                   #
################################################

########## imshow [im1, im2] ###################

def imagesc(*args):

    N = len(args)

    if N == 0 or N > 6:
        print ("pyraft:: imagesc error: max of 4 images")
        return
    
    if N == 1:
        imshow(args[0], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        show()
        return
    
    if N == 2:
        plt.subplot(1, 2, 1)
        imshow(args[0], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(1, 2, 2)
        imshow(args[1], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        show()
        return
        
    if N == 3:
        plt.subplot(1, 3, 1)
        imshow(args[0], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(1, 3, 2)
        imshow(args[1], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(1, 3, 3)
        imshow(args[2], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        show()
        
        return
        
    if N == 4:
        plt.subplot(2, 2, 1)
        imshow(args[0], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 2, 2)
        imshow(args[1], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 2, 3)
        imshow(args[2], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 2, 4)
        imshow(args[3], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        show()
        
        return       

    if N == 5:
        plt.subplot(2, 3, 1)
        imshow(args[0], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 3, 2)
        imshow(args[1], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 3, 3)
        imshow(args[2], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 3, 4)
        imshow(args[3], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 3, 5)
        imshow(args[4], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        show()
        
        return       

    if N == 6:
        plt.subplot(2, 3, 1)
        imshow(args[0], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 3, 2)
        imshow(args[1], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 3, 3)
        imshow(args[2], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 3, 4)
        imshow(args[3], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 3, 5)
        imshow(args[4], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        plt.subplot(2, 3, 6)
        imshow(args[5], extent=[0, 100, 0, 1], aspect=100)
        colorbar()
        show()
        
        return       



def plotsc(y):
    plt.figure()
    plt.plot(y)
    plt.show()
    return 


    
