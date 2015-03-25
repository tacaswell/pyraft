import numpy
from scipy.ndimage.interpolation import rotate, shift
from ..raftypes import image

######Functions for camera simulation in snapshot
def snapshot_affine(img, row):	
	img_rot = numpy.zeros(img.shape)
	img_rot = rotate(img, row, reshape=0)
	return img_rot

def snapshot_shift(img, dx, dy):	
	img_shift = numpy.zeros(img.shape)
	img_shift = shift(img, [dy, dx])
	return img_shift

######Functions for single slice sample simulation
def image_rotate(img, row):
	img_rot = numpy.zeros(img.shape)
	img_rot = rotate(img, row, reshape=0)
	img_rot = image(img_rot)
	img_rot.bottom_right = img.bottom_right
	img_rot.top_left = img.top_left
	return img_rot

def image_shift(img, dx, dy):
	img_shift = numpy.zeros(img.shape)
	img_shift = shift(img, [dy, dx])
	img_shift = image(img_shift)
	img_shift.bottom_right = img.bottom_right
	img_shift.top_left = img.top_left
	return img_shift


######Functions for block sample simulation
def vol_rotate(block, theta, row, pitch, pivot):
	sizeX = block.shape[0]
	sizeY = block.shape[1]
	sizeZ = block.shape[2]
	
	block_new = numpy.zeros(block.shape)	

	for cut in numpy.arange(0, sizeY):
		block_new[:, cut, :] = rotate(block[:, cut, :], theta, reshape=0)

	for cut in numpy.arange(0, sizeX):
		block_new[cut, :, :] = rotate(block_new[cut, :, :], row, reshape=0)

	for cut in numpy.arange(0, sizeZ):
		block_new[:, :, cut] = rotate(block_new[:, :, cut], pitch, reshape=0)

	return block_new


def vol_shift(block, dx, dy, dz):
	block_new = numpy.zeros(block.shape)	
	block_new = shift(block, [dy, dx, dz])
	return block_new



def rotateImage(img, angle, pivot):
    padX = [img.shape[1] - pivot[0], pivot[0]]
    padY = [img.shape[0] - pivot[1], pivot[1]]
    imgP = numpy.pad(img, [padY, padX], 'constant')
    imgR = rotate(imgP, angle, reshape=False)
    return imgR[padY[0] :-padY[1], padX[0] :-padX[1]]
