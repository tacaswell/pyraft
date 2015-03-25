import numpy

###################################################################
# #Import pyraft with full wrappers
import pyraft

###################################################################
# #Generate simple phantom, sinogram, backprojection and filtered backprojection
f = pyraft.gear()
s = pyraft.radon(f)
b = pyraft.backprojection(s)
x = pyraft.fbp(s)

# notice the standard shape sizes
print f.shape, s.shape, b.shape, x.shape

# Show the images
pyraft.imagesc(f, s, b, x)

###################################################################
# #Generate a simple volume, a snapshot and theta=30 snapshot
v = pyraft.vol_gear()

# Create the snapshots
snap = pyraft.snapshot(v)
snap30 = pyraft.snapshot(v, 30 * numpy.pi / 180)

# Show the images
pyraft.imagesc(snap, snap30)


###################################################################
# #Compares Snapshot of 3D phantom vs 2D phantom
# it is important to notice that the phantom gear is exactly the snapshot at theta=0 (both normalized)
v = pyraft.vol_gear()
snap = pyraft.snapshot(v)
snap = snap / snap.max().max()
f = pyraft.gear([128, 128])
f = f / f.max().max()

# Show the images
pyraft.imagesc(snap, f, snap - f)

###################################################################
# #Test pyraft center sino function (pre-Radon)


center = 10
x = pyraft.gear()
s = pyraft.radon(x)
s2 = pyraft.image_shift(s, 0 , center)
noise = abs(numpy.random.poisson(.05, s2.shape))
s3 = s2 + noise

s_center, offset = pyraft.center_sino(s2)
error = (center - offset) / center
print offset
print 'Error of ' + str(abs(error * 100)) + '%'

i = pyraft.fbp(s)
i2 = pyraft.fbp(s2)
i_center = pyraft.fbp(s_center)

pyraft.imagesc(s, s2, s_center, i, i2, i_center)
pyraft.imagesc(s - s2, s - s_center, i - i2, i - i_center)



s_center, offset = pyraft.center_sino(s3)
error = (center - offset) / center
print offset
print 'Noise Error of ' + str(abs(error * 100)) + '%'

i = pyraft.fbp(s)
i2 = pyraft.fbp(s3)
i_center = pyraft.fbp(s_center)

pyraft.imagesc(s, s3, s_center, i, i2, i_center)
pyraft.imagesc(s - s3, s - s_center, i - i2, i - i_center)



###################################################################
# #Test pyraft EM cap.

f = pyraft.gear()

s = pyraft.radon(f, [367, 180])
noise = numpy.random.poisson(1, s.shape)
s1 = s + noise


a10 = pyraft.em(s1, [256, 256], niter=10)

pyraft.imagesc(f, a10, f - a10)


###################################################################
# #Test pyraft GPU cap.

f = pyraft.gear([2048, 2048])
s = pyraft.radon_gpu(f, 2048, 200)
b = pyraft.back_gpu(s)
x = pyraft.fbp_gpu(s)
# notice the standard shape sizes
print f.shape, s.shape, b.shape, x.shape

# Show the images
pyraft.imagesc(f, s, b, x)

xem = pyraft.em_gpu(s,2048,20)

pyraft.imagesc(f,x,xem)
