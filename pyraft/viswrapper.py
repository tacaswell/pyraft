from mayavi import mlab


#############

def volRender(vol):
	mlab.figure(bgcolor=(1, 1, 1))
	mlab.pipeline.volume(mlab.pipeline.scalar_field(vol), vmin=0, vmax=vol.max().max())
	mlab.outline()
	mlab.show()	

def volContour(vol):
	mlab.figure(bgcolor=(1, 1, 1))
	

	mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(vol),
                                 plane_orientation='x_axes',
                                 slice_index=10,
                                 )

	mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(vol),
                                 plane_orientation='y_axes',
                                 slice_index=10,
	                         )

	mlab.outline()
	mlab.show()

