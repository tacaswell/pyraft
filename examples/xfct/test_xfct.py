import numpy
import pyraft

dens = numpy.transpose( pyraft.image( numpy.loadtxt('new.dat') ) )
trans = numpy.transpose( pyraft.image( numpy.loadtxt('ein.dat') ) )
fluor = numpy.transpose( pyraft.image( numpy.loadtxt('hilb.dat') ) )

pyraft.imagesc(dens,trans,fluor)

d = pyraft.xfct.radon(dens, trans, fluor, [dens.shape[0], 180])

b = pyraft.xfct.backprojection360(d)

b2 = pyraft.xfct.backprojection(d, trans, fluor, dens.shape)

rec = pyraft.xfct.akt(d, trans, fluor, dens.shape, 10, 0)
#rec = pyraft.xfct.em(d, trans, fluor, dens.shape, 50, 0)
#pyraft.imagesc(dens, rec)

y = pyraft.xfct.fbp360( d )

pyraft.imagesc(b,b2,y, rec)
