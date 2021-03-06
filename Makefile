RUNDIR=	raft xfct pyraft pyraft/sinogram pyraft/phantom pyraft/backprojection pyraft/filters pyraft/inversion pyraft/iterative pyraft/coordinates pyraft/xfct

all: build install

build: 
	sudo python setup.py build

install:
	sudo -E python setup.py install

clean:
	sudo python setup.py clean --all
	rm -fr build/ *.egg-info/ dist/	*~
	@for j in ${RUNDIR}; do rm -rf $$j/*.pyc; rm -rf $$j/*~; done

