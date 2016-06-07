###########################################################################
#
#  CENTERING SINOGRAMS: FUNCTIONS
#
###########################################################################

__author__     = "Eduardo X Miqueles"
__copyright__  = "Copyright 2015, CNPEM/LNLS" 
__credits__    = ""
__maintainer__ = "Eduardo X Miqueles"
__email__      = "edu.miqueles@gmail.com" 
__date__       = "23.Aug.2015" 


import sys
from scipy import ndimage
from scipy.misc import imsave
from scipy import signal as scipy_signal
from scipy.optimize import minimize as scipy_minimize
from scipy.optimize import leastsq
from scipy.special import comb as nchoosek

import matplotlib.pyplot as plt
import numpy

###########################################################################
#
# Overview:	Azevedo's analytical inverse matrix 3x3
#		
#		Only upper part of the matrix is needed
#
#		TODO:	complete lower part using 3 cofactor matrices
#			in case we need to compute center os mass
#	
#		DOI: 10.1109/23.55866	 
###########################################################################

def azevedo_matrix(n):

    H = numpy.zeros([3,3])

    n24 = (float(n*n)/4)
    n22 = (float(n*n)/2)
    n2  = (float(n)/2)

    pi2n = (numpy.pi/(2*n))
    t = numpy.tan( pi2n )
    t1 = (1.0/t)
    d = n2*(n22 - 1 - t1**2)
    d1  = (1.0/d)

    H[0,0] = d1*n24
    H[0,1] = -d1*n2
    H[0,2] = -d1*n2*t1
    H[1,1] = d1*( n22 - t1**2)
    H[1,2] = d1*t1
    H[2,2] = d1*(n22 - 1)
    H[1,0] = H[0,1]
    H[2,0] = H[0,2]
    H[2,1] = H[1,2]
    
    return H

###########################################################################
#
# Overview:	Azevedo's analytical solution for shift
#		
#		DOI: 10.1109/23.55866	 
###########################################################################

def azevedo_analytical_shift(t0, tf, sino):
   
    N = sino.shape[1]
    R = sino.shape[0]
 
    th = numpy.linspace(0, 180, N, endpoint=False) * (numpy.pi/180)
    th.shape = [len(th), 1]
    t = numpy.linspace(t0, tf, R)
    t.shape = [len(t), 1]
    t0 = float(t[0])   
    dt = float(t[1] - t0)
 
    a1 = numpy.ones([N, 1])
    a2 = numpy.cos(th)
    a3 = numpy.sin(th)
 
    A = numpy.hstack([a1, a2, a3])
 
    T = numpy.kron(t, numpy.ones([1, N]))
    T = numpy.flipud(T) #because t is flipped

    m = sino.sum(0)
    mask = numpy.abs(m) < 1e-5    
    m[mask] = 1.0
 
    s1 = T * sino
    s = s1.sum(0) / m
    s.shape = [len(s), 1]
 
    H = azevedo_matrix(N)
  
    vec = numpy.dot( H, numpy.dot(numpy.transpose( A ), s) )
    
    return -float(vec[0]), float(vec[1]), float(vec[2])  
    ## -beta, x_cmass, y_cmass
    ## -beta because of my definition of shift on t-axis    

###########################################################################
#
# Overview:	Find the offset of a given 'value' at a mesh in [t0,tf] with
#		N points. Compute the offset for each double in the input
#		array 'value'
#
###########################################################################

def get_offset_array(value, t0, tf, N):
     
    t = numpy.linspace(t0, tf, N)
    t.shape = [len(t), 1]
    dt = float(t[1] - t0)
 
    ind_zero = ( numpy.ceil( (0-t0)/dt ) )
    ind_value =  ( numpy.ceil( (value-t0)/dt ) )
 
    offset = ind_zero - ind_value
 
    return offset
 
###########################################################################
#
# Overview:	Exactly the same as function 'get_offset_array', but
#		input is just a double, not an array
#		
###########################################################################

def get_offset(value, t0, tf, N):
     
    t = numpy.linspace(t0, tf, N)
    t.shape = [len(t), 1]
    dt = float(t[1] - t0)
 
    ind_zero = ( numpy.ceil( (0-t0)/dt ) )
    ind_value =  ( numpy.ceil( (value-t0)/dt ) )
 
    offset = int(ind_zero - ind_value)
 
    return offset

###########################################################################
#
# Overview:	Translate a signal y(t) to y(t-c) using FFT
#		
###########################################################################
 
def translate(y, c, t0, tf):
    #translate array y(t) to y(t-c)
 
    dt = (tf-t0)/(float(len(y)))
    wc = 1.0/(2.0*dt)
    w = numpy.linspace(-wc, wc, len(y))
    sig = numpy.exp(-2* numpy.pi * 1j * c * w)   
  
    fy = numpy.fft.fft(y)
    w = numpy.fft.ifft( fy * numpy.fft.fftshift(sig) ) 
 
    return w.real

###########################################################################
#
# Overview:	Zero padding of a sinogram
#		
#		TODO: remove this function in the future
###########################################################################
 
def zpSino(sino):

    N = sino.shape[1]
    R = sino.shape[0]  
 
    ssino = numpy.zeros([200+R,N]) 
    ssino[0:100,:] = numpy.zeros([100,N])
    ssino[100:100+R,:] = sino
    ssino[100+R:200+R,:] = numpy.zeros([100,N])
     
    return ssino

###########################################################################
#
# Overview:	Equivalent to numpy.roll() function, but using subpixel 
#		precision through FFT. Design only for sinogram images.
#		(rays x angles)
#
#		TODO: make it fast
###########################################################################
 
def rollSino(sino, c):
 
	N = sino.shape[1]

	new = numpy.zeros(sino.shape)
	for j in range(N): 
		y = numpy.flipud(sino[:,j])
		w = translate(y, c, -1.0, 1.0) 
		new[:,j] = numpy.flipud(w)

	return new
 
###########################################################################
#
# Overview:	Functions for the problem: x = argmin[ func4; params ]
#
###########################################################################

#
# auxiliary function: compute a given matrix A[k]
#

def matrix_A(k, th):
	V = len(th)
	if k == 0:
		A = numpy.ones((V,1))
	else:
		c = numpy.cos(th)
		s = numpy.sin(th)
		A = numpy.zeros((V,k+1))
	for j in range(0,k+1):
		vec = (c**(k-j)) * (s**(j)) * nchoosek(k,j)
		vec.shape = [V,]
		A[:,j] = vec
	return A

#
# auxiliary function: compute the constant vector b[k]
#

def func4_data(sino, k, t0, tf):
  
    V = sino.shape[1]
    R = sino.shape[0]
 
    th = numpy.linspace(0, 180, V, endpoint=False) * (numpy.pi/180)
    th.shape = [len(th), 1]
    t = numpy.linspace(t0, tf, R)
    t.shape = [len(t), 1]
     
    
    m = sino.sum(0)
    mask = numpy.abs(m) < 1e-5    
    m[mask] = 1.0

    T = numpy.kron(t**k, numpy.ones([1, V]))
    T = numpy.flipud(T) #because t is flipped
 
    s = T * sino
   
    b = s.sum(0) / m
    b.shape = [len(b), 1]
     
    return b

#
# auxiliary function: compute a the diagonal matrix D(beta)
#
	
def func4_matrix_diag(beta, k, V):
       
    diag_ = beta**k
    diff_diag_ = k*beta**(k-1)
 
    for j in range(1,k+1):
	    temp_diag = numpy.ones((j+1,)) * nchoosek(k,k-j) * (beta**(k-j))
	    diag_ = numpy.r_[ temp_diag, diag_ ] 
	    if (k-j) == 0:
		    temp_diag = numpy.zeros((j+1,))
		    diff_diag_ = numpy.r_[ temp_diag, diff_diag_ ]
	    else:
		    temp_diag = numpy.ones((j+1,)) * nchoosek(k,k-j) * (k-j) * (beta**(k-j-1))
		    diff_diag_ = numpy.r_[ temp_diag, diff_diag_ ]
 	
    D  = numpy.diag(diag_)
    dD = numpy.diag(diff_diag_)  

    dim = len(diag_)
    U = numpy.zeros([dim-1, dim-1])
    U[:,:] = D[0:dim-1,0:dim-1]

    dU = numpy.zeros([dim-1, dim-1])
    dU[:,:] = dD[0:dim-1,0:dim-1]

    return D, U, dD, dU

#
# auxiliary function: compute block matrix A
#

def func4_matrix_blocks(k, t0, tf, V):
    th = numpy.linspace(0, 180, V, endpoint=False) * (numpy.pi/180)
    th.shape = [len(th), 1]
    A = matrix_A(0, th)
    for j in range(1,k+1):
        temp_A = matrix_A(j, th)
        A = numpy.c_[temp_A, A]
    dim = A.shape[1]
    S = numpy.zeros([V, dim-1])
    S[:,:] = A[:,0:dim-1]
    return A, S

#
# function to minimize
#

def func4(x, *args):
    
    # -- remember: x = (beta, y)

    params = args[0]
    sino = params[0]
    k = params[4] 
    A = params[5]
    b = params[6]
     
    V = sino.shape[0]
    
    beta = x[0] 
      
    # -- diagonal depending on beta
    
    D, _, _, _ = func4_matrix_diag(beta, k, V)
    
    # -- 
    
    M = numpy.dot( A, D )

    y = x[1:len(x)]
    y.shape = [len(y), 1]
    z = numpy.r_[ y,  numpy.ones([1,1]) ]

    w = b - numpy.dot(M, z)
    
    fun = 0.5 * (numpy.linalg.norm(w))**2
  
    return fun

#
# gradient of function: func4()
#

def jac_func4(x, *args):
    
    # -- remember: x = (beta, y)

    params = args[0]

    sino = params[0]
    k = params[4] 
    A = params[5]
    b = params[6]
    S = params[7]
     
    V = sino.shape[1]
    
    beta = x[0] 

    D, U, dD, dU = func4_matrix_diag(beta, k, V)

    y = x[1:len(x)]
    y.shape = [len(y), 1]
    
    e = numpy.ones([V,1])
    L = numpy.dot(S, U)
    LtLy = numpy.dot( numpy.dot(numpy.transpose(L), L), y)
    Lte = numpy.dot(numpy.transpose(L), e)
    
    # diff: variable y
    grad_y = LtLy + (beta**k) * Lte
 
    # diff: variable beta

    M = numpy.dot( A, D )

    z = numpy.r_[ y, numpy.ones([1,1]) ]

    dMz = numpy.dot( numpy.dot(A,dD), z) 
    grad_beta =  - numpy.dot( numpy.transpose(b - numpy.dot(M,z)), dMz)
 
    return numpy.r_[ grad_y, grad_beta ]


def func5(x, *args):
    
    # -- remember: x = (beta, y)

    params = args[0]
    sino = params[0]
    k = params[4] 
    M = params[5]
    b = params[6]
    y = params[7]    
    beta = params[8] 
 
    z = numpy.zeros([k+1 + len(y) + 1,])
    z[0:k+1] = x[:]
    z[k+1:k+1+len(y)] = y[:,0]
    z[k+1 + len(y)] = 1
    z.shape = [len(z), 1]

    w = b - numpy.dot(M, z)

    fun = 0.5 * (numpy.linalg.norm(w))**2
  
    return fun


###########################################################################
#
# Overview:	Consistency function
#		
#		Given a sinogram S (rays x angles) and a shift x, compute 
#		the relative difference - using an L2 norm:
#		
#		rest =  | shift(a,x) - shift(b,x) | / |a|,
#  
#		where a = S[:,pi] and b = flipud( S[:,0])
# 
###########################################################################

def consistency(x, *args):
    # input sinogram: (rays x angles)
 
    params = args[0]
 
    sino = params[0]
    epsilon = params[1]  
    t0 = params[2]
    tf = params[3]   
   
    #print sino.shape, epsilon, t0, tf

    N = sino.shape[1]
    R = sino.shape[0]
 
    a = sino[:,N-1]
    b = numpy.flipud(sino[:,0])
    y = translate(a, -x[0], t0, tf) - translate(b, x[0], t0, tf)
         
    rest = (numpy.linalg.norm(y))/(numpy.linalg.norm(a))
    return rest**2

    #rest = - (numpy.linalg.norm(y)/numpy.linalg.norm(a)) + epsilon
    #return rest


###########################################################################
#
# Overview:	Phase-Correlation
#		
#		Given a sinogram S (rays x angles) and a shift x, compute 
#		the positive phase correlation using FFT between a and b
#		where a = S[:,pi] and b = flipud( S[:,0])
# 
###########################################################################

def correlation(sino):
    # input sinogram: (rays x angles)
  
    N = sino.shape[1]
    R = sino.shape[0]
 
    a = sino[:,N-1]
    b = numpy.flipud(sino[:,0])
  
    t = numpy.linspace(-1., 1., R)    

    xcorr = numpy.correlate(a,b,"same")

    idx = numpy.argmax(xcorr)
    
    beta = t[idx]

    return beta

###########################################################################
#
# Overview:	Consistency constraint
#		
#		Given a sinogram S (rays x angles) and a shift x, compute 
#		the positive constraint rest(x) >= 0 - using an L2 norm:
#		
#		rest =  - | shift(a,x) - shift(b,x) | + epsilon,
#  
#		where a = S[:,pi] and b = flipud( S[:,0])
# 
###########################################################################

def consistency_cstr(x, *args):
    # input sinogram: (rays x angles)
 
    params = args[0]
 
    sino = params[0]
    epsilon = params[1]  
    t0 = params[2]
    tf = params[3]
 
    N = sino.shape[1]
    R = sino.shape[0]
 
    a = sino[:,N-1]
    b = numpy.flipud(sino[:,0])
    y = translate(a, x[0], t0, tf) - translate(b, -x[0], t0, tf)

    #rest =  - (numpy.linalg.norm(y)/numpy.linalg.norm(a))**2 + epsilon	
    
    rest =  - (numpy.linalg.norm(y)/numpy.linalg.norm(a)) + epsilon

    return rest**2

###########################################################################
#
# OFFSET Algorithm #1
#
# Overview:	Azevedo's analytical algorithm
#		DOI: 10.1109/23.55866 
#
###########################################################################

def offset1(sino):
 
    t0 = -1.0
    tf =  1.0

    beta, x, y = azevedo_analytical_shift(t0, tf, sino)
 
    offset = get_offset( beta, t0, tf, sino.shape[0])

    return beta, offset, numpy.array([x,y])

####################################################################
#
# OFFSET Algorithm #2
#
# Overview:	Prince's least square algorithm
#		DOI: 10.1109/83.236529
#
#
####################################################################

def offset2(sino):
 
    N = sino.shape[1]
    R = sino.shape[0]
 
    th = numpy.linspace(0, 180, N, endpoint=False) * (numpy.pi/180)
    th.shape = [len(th), 1]
    t = numpy.linspace(-1.0, 1.0, R)
    t.shape = [len(t), 1]
 
    a2 = numpy.cos(th)
    a3 = numpy.sin(th)
 
    A = numpy.hstack([a2, a3])
 
    T = numpy.kron(t, numpy.ones([1, N]))
    T = numpy.flipud(T) #because t is flipped 

    m = sino.sum(0)
    mask = numpy.abs(m) < 1e-5
    m[mask] = 1.0
   
    w = T * sino
    b = w.sum(0)/ m
    b.shape = [N, 1]
 
    z = numpy.linalg.lstsq(A, b)[0]
     
    x = float(z[0])
    y = float(z[1])
 
    fit2 = numpy.dot(A,z)
    pixel = b - fit2
    pixel = - pixel  ## shift because my definition of rolling on t-axis

    offset_arr = get_offset_array(pixel, -1.0, 1.0, R) 
     
    d = N-1
    offset = int( offset_arr[d] )
    c = float(pixel[d])
 
    newsino = numpy.zeros( sino.shape )
    for j in range(0,sino.shape[1]): 
        newsino[:,j] = numpy.roll( sino[:,j], int(offset_arr[j]) )

    return c, offset, newsino


####################################################################################
#
# OFFSET Algorithm #3
#
# Overview:	Minimization of consistency (0-180 degrees correlation). 
#		=> Brute force code :-(
#
###################################################################################

def offset3( sino, raymax, start ):

    epsilon = 1e-8
    t0 = -1.0
    tf =  1.0
     
    params = (sino, epsilon, t0, tf)
 
    z0 = [start]
    
    cons = ({'type': 'ineq', 'fun' : lambda x: - x[0] + raymax },
            {'type': 'ineq', 'fun' : lambda x:   x[0] + raymax })


    ###

    sol_con = scipy_minimize(consistency, z0, args=(params,), constraints=cons)
    beta = sol_con.x[0]

    ###

    offset = get_offset( beta, -1.0, 1.0, sino.shape[0])

    return beta, offset, offset
    

####################################################################
#
# OFFSET Algorithm #4 (new)
#
# Overview:	Solve the problem x = argmin [ func4; params]
#
#		s.t.: 		consistency_cstr(x) >= 0
#				| x | <= RAYMAX 		
#
####################################################################

def offset4(sino, k):
	
    V = sino.shape[1]
    R = sino.shape[0]

    MAX_NPIXELS = float(R/4)

    epsilon = 1e-8
    t0 = -1.0
    tf =  1.0

    dt = (tf-t0)/R

    RAYMAX = MAX_NPIXELS * dt

    dimension = int( (int(k)+1)*(int(k)+2)/2 )

    b = func4_data(sino, k, t0, tf)
    A,S = func4_matrix_blocks(k, t0, tf, V)

    z0 = numpy.zeros((dimension,1)) 
    z0[0,0], z0[1,0], z0[2,0] = azevedo_analytical_shift(t0, tf, sino)  

    # -- optimization using scipy

    params = (sino, epsilon, t0, tf, k, A, b, S)

    cons = ({'type': 'ineq', 'fun' : lambda x, p: consistency_cstr(x, p), 'args': (params,)},
            {'type': 'ineq', 'fun' : lambda x: - x[0] + RAYMAX },
            {'type': 'ineq', 'fun' : lambda x: x[0] + RAYMAX })
 
    solution = scipy_minimize(func4, z0, args=(params,), constraints=cons, method='SLSQP')

    ## --- 

    beta = solution.x[0]

    offset = get_offset( beta, -1.0, 1.0, sino.shape[0])

    return beta, offset, solution.x[1:dimension]


####################################################################
#
# OFFSET Algorithm #5 (new)
#
# Overview:	Solve the problem x = argmin [ func4; params]
#
####################################################################

def offset5(sino, k):
	    
    V = sino.shape[1]
    R = sino.shape[0]

    MAX_NPIXELS = 400

    epsilon = 1e-5
    t0 = -1.0
    tf =  1.0

    dt = (tf-t0)/R

    RAYMAX = MAX_NPIXELS * dt

    dimension = int( (int(k)+1)*(int(k)+2)/2 )

    b = func4_data(sino, k, t0, tf)
    A,S = func4_matrix_blocks(k, t0, tf, V)

    z0 = numpy.zeros((dimension,1))
    z0[0,0], _, _ = azevedo_analytical_shift(t0, tf, sino)

    # -- optimization using scipy

    params = (sino, epsilon, t0, tf, k, A, b, S)
    
    sol_con = scipy_minimize(func4, z0, args=(params,))

    ## --- 

    beta = sol_con.x[0]

    offset = get_offset( beta, -1.0, 1.0, sino.shape[0])

    return beta, offset, sol_con.x[1:dimension]


####################################################################
#
# OFFSET Algorithm #6 (new)
#
# Overview:	Solve the problem x = argmin [ func4; params] 
#		recursively
#
# Prototype:	offset( sino, k, (percentage_rays, initial_point) )
####################################################################

def offset6(sino, k, *args):

    #
    # parameters to speed up method
    # 
    if len(args)>0:
    	argum = args[0]
    ############    

    V = sino.shape[1]
    R = sino.shape[0]

    th = numpy.linspace(0, 180, V, endpoint=False) * (numpy.pi/180)
    th.shape = [len(th), 1]

    if len(args) and len(argum)==1:
	#MAX_NPIXELS = int(argum[0]*R)
        RAYMAX = float(argum[0])
    else:
	#MAX_NPIXELS = int(0.25*R)
        RAYMAX = 0.5

    epsilon = 1e-8
    t0 = -1.0
    tf =  1.0

    dt = (tf-t0)/R

    #RAYMAX = MAX_NPIXELS * dt

    moments = numpy.zeros([k+1, k+1])
    moments[0][0] = 1.0
	
    # -- first optimization step

    j = 1

    b = func4_data(sino, j, t0, tf)
    A,S = func4_matrix_blocks(j, t0, tf, V)

    bnds = ((-RAYMAX, RAYMAX),(None,None),(None,None))

    params = (sino, epsilon, t0, tf, j, A, b, S)

    cons = ({'type': 'ineq', 'fun' : lambda x, p: consistency_cstr(x, p), 'args': (params,)},
            {'type': 'ineq', 'fun' : lambda x: - x[0] + RAYMAX },
            {'type': 'ineq', 'fun' : lambda x: x[0] + RAYMAX })

    dimension = int( (int(j)+1)*(int(j)+2)/2 )

    z0 = numpy.zeros((dimension,1))
    beta0, cxmass, cymass = azevedo_analytical_shift(t0, tf, sino)
        

    if len(args)>0 and len(argum)==2:
        initial = argum[1]
        z0[0,0] = initial
        z0[1,0] = cxmass
        z0[2,0] = cymass
    else:
        z0[0,0] = -beta0
        z0[1,0] = cxmass
        z0[2,0] = cymass

    solution = scipy_minimize(func4, z0, args=(params,), constraints=cons)#, bounds=bnds) # method='SLSQP')

    beta = solution.x[0]

    y = solution.x[1:dimension]
    y.shape = [len(y), 1]
    
    moments[1,0] = y[0] # z0[1,0] #y[0]
    moments[0,1] = y[1] # z0[2,0] #y[1]

    if k > 1:

    	for m in range(2,k+1):
            dimension = m+1
            b = func4_data(sino, m, t0, tf)
            D, _, _, _ = func4_matrix_diag(beta, m, V)
            A_ = matrix_A(m, th)
            A = numpy.c_[ A_, A ]
            M = numpy.dot(A, D)	
            z0 = numpy.zeros((dimension,1)) 

            # -- optimization using scipy

            params = (sino, epsilon, t0, tf, m, M, b, y, beta)
            
            solution = scipy_minimize(func5, z0, args=(params,))

            y_ = solution.x
            
            y_.shape = [dimension, 1]

            y = numpy.r_[ y_, y]

            ## --- moment matrix
            
            for l in range(0,m+1):
                    moments[m-l][l] = y[l]
    

    offset = get_offset( beta, -1.0, 1.0, sino.shape[0])
           
    return -beta, offset, moments


##############################################################
#
# MAIN FUNCTION(S) FOR PYRAFT
#
##############################################################

def centerSino( sino, k, *args ):

    #input sinogram: (rays x nangles)

    return offset6( sino, k, *args )


def centering( sino, maxoffset, *args):

    #input sinogram: (rays x nangles)

    if len(args)>0:
	#comparison between input values: given by user and correlation
	#choice: the one giving lower value for the objective function.
        argum = args[0]
        shift_ = numpy.zeros([2,])
        retorno = numpy.zeros([2,])
        fun = numpy.zeros([2,])

        starting = numpy.array([argum, correlation(sino)])      

        retorno[0],shift_[0],_ = offset3( sino, maxoffset, starting[0])
        retorno[1],shift_[1],_ = offset3( sino, maxoffset, starting[1])

        fun[0] = consistency( [retorno[0]], (sino, 1e-8, -1.,1.))
        fun[1] = consistency( [retorno[1]], (sino, 1e-8, -1.,1.))	

        idx = numpy.argmin(fun)

        beta = retorno[idx]         
        shift = shift_[idx]
   
    else: 

        starting = correlation(sino)
 
        beta,shift,_ = offset3( sino, maxoffset, starting) 
 
    ###############
	
    return beta, shift, shift


