#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from os.path import join as pjoin
import warnings
import glob
from setuptools import setup, Extension, find_packages
from distutils.extension import Extension
from distutils.command.build_ext import build_ext
import subprocess
import numpy
import sys

##########################################################

compile_xfct=0

if '--xfct' in sys.argv:
    compile_xfct = 1
    print('Compiling XFCT functions with GSL & Confuse')
    sys.argv.remove('--xfct')
    
# Set Python package requirements for installation.
install_requires = [
            'numpy>=1.8.0',
            'scipy>=0.14.0',
            'h5py>=2.0',
	    'pynfft>=1.3.2'
            ]


def find_in_path(name, path):
    "Find a file in a search path"
    #adapted fom http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None


def locate_cuda():
    """Locate the CUDA environment on the system

    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDAHOME env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """

    # first check if the CUDAHOME env variable is in use
    if 'CUDAHOME' in os.environ:
        home = os.environ['CUDAHOME']
        nvcc = pjoin(home, 'bin', 'nvcc')
    else:
        # otherwise, search the PATH for NVCC
        nvcc = find_in_path('nvcc', os.environ['PATH'])
        if nvcc is None:
            print ('The nvcc binary could not be '
                'located in your $PATH. Either add it to your path, or set $CUDAHOME')
            return None
        home = os.path.dirname(os.path.dirname(nvcc))

   
    check = pjoin(home, 'lib')

    if not os.path.exists(check):
	cudaconfig = {'home':home, 'nvcc':nvcc,
        	      'include': pjoin(home, 'include'),
                      'lib': pjoin(home, 'lib64')}
    else:
	cudaconfig = {'home':home, 'nvcc':nvcc,
        	      'include': pjoin(home, 'include'),
                      'lib': pjoin(home, 'lib')}


    for k, v in cudaconfig.iteritems():
        if not os.path.exists(v):
            print ( 'The CUDA %s path could not be located in %s' % (k, v))
    return cudaconfig

CUDA = locate_cuda()

# enforce these same requirements at packaging time
import pkg_resources
for requirement in install_requires:
    try:
        pkg_resources.require(requirement)
    except pkg_resources.DistributionNotFound:
        msg = 'Python package requirement not satisfied: ' + requirement
        msg += '\nsuggest using this command:'
        msg += '\n\tpip install -U ' + requirement.split('=')[0].rstrip('>')
        print (msg)
        raise (pkg_resources.DistributionNotFound)

#########################################
# set XFCT libraries @ pyraft/ratypes.py

import shutil

if compile_xfct:
    with open("pyraft/raftypes.py", "wt") as fout:
        with open("pyraft/raftypes_basis.py", "rt") as fin:
            for line in fin:
                fout.write(line.replace('#XFCT_LIBS', 'libgsl=ctypes.CDLL( ctypes.util.find_library( "lgsl" ), mode=ctypes.RTLD_GLOBAL )'+str('\n')+'libgslcblas = ctypes.CDLL( ctypes.util.find_library( "lgslcblas" ), mode=ctypes.RTLD_GLOBAL )'+str('\n')+'libconfuse  = ctypes.CDLL( ctypes.util.find_library( "lconfuse" ), mode=ctypes.RTLD_GLOBAL)'))
else:
    shutil.copyfile("pyraft/raftypes_basis.py", "pyraft/raftypes.py")


########################################################

if CUDA:
	if compile_xfct:
		raft_codes = list( set(glob.glob('raft/raft_*.c*')) | set(glob.glob('xfct/*.c')) )	 
	else:
		raft_codes = set(glob.glob('raft/raft_*.c*'))

else:
	if compile_xfct:
		raft_codes = list( (set(glob.glob('raft/raft_*.c*'))-set(glob.glob('raft/*.cu'))) | set(glob.glob('xfct/*.c')) )
	else:
		raft_codes = set(glob.glob('raft/raft_*.c*'))-set(glob.glob('raft/*.cu'))


# Create reconstruction shared-library.
if CUDA:
	if compile_xfct:
            ext_raft = Extension(name='pyraft.lib.libraft',
                                 sources=list(raft_codes),
                                 library_dirs=[CUDA['lib']],
                                 libraries=['cudart'],
                                 runtime_library_dirs=[CUDA['lib']],
                                 extra_compile_args={'gcc': ['-pedantic','-std=c++0x'],
                                                     'nvcc': ['-Xcompiler','-use_fast_math', '--ptxas-options=-v', '-c', '--compiler-options', '-fPIC']},
                                 extra_link_args=['-std=c++0x','-lfftw3_threads','-lfftw3','-lm','-lblas','-lpthread','-lgsl','-lgslcblas','-lconfuse'],		
                                 include_dirs = [ CUDA['include'], ])
        else:
            ext_raft = Extension(name='pyraft.lib.libraft',
                                 sources=list(raft_codes),
                                 library_dirs=[CUDA['lib']],
                                 libraries=['cudart'],
                                 runtime_library_dirs=[CUDA['lib']],
                                 extra_compile_args={'gcc': ['-pedantic','-std=c++0x'],
                                                     'nvcc': ['-Xcompiler','-use_fast_math', '--ptxas-options=-v', '-c', '--compiler-options', '-fPIC']},
                                 extra_link_args=['-std=c++0x','-lfftw3_threads','-lfftw3','-lm','-lblas','-lpthread'],		
                                 include_dirs = [ CUDA['include'], ])
else:
    if compile_xfct:
        ext_raft = Extension(name='pyraft.lib.libraft',
		             sources=list(raft_codes),
		             extra_compile_args={'gcc': ['-pedantic','-std=c++0x']},
		             extra_link_args=['-std=c++0x','-lfftw3_threads','-lfftw3','-lm','-lblas','-lpthread','-lgsl','-lgslcblas','-lconfuse'])
    else:
        ext_raft = Extension(name='pyraft.lib.libraft',
		             sources=list(raft_codes),
		             extra_compile_args={'gcc': ['-pedantic','-std=c++0x']},
		             extra_link_args=['-std=c++0x','-lfftw3_threads','-lfftw3','-lm','-lblas','-lpthread'])



def customize_compiler_for_nvcc(self):
    """inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.
    
    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on."""
    
    # tell the compiler it can processes .cu
    self.src_extensions.append('.cu')

    # save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            # use the cuda for .cu files
            self.set_executable('compiler_so', CUDA['nvcc'])
            # use only a subset of the extra_postargs, which are 1-1 translated
            # from the extra_compile_args in the Extension class
            postargs = extra_postargs['nvcc']
        else:
            postargs = extra_postargs['gcc']

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # inject our redefined _compile method into the class
    self._compile = _compile


# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)


# Main setup configuration.
setup(
    name='pyraft',
    version=open('VERSION').read().strip(),
    
    packages = find_packages(),
    include_package_data = True,
    
    ext_modules=[ext_raft],
      # inject our custom trigger
    cmdclass={'build_ext': custom_build_ext},

      # since the package has c code, the egg cannot be zipped
    zip_safe=False,    

    author='Eduardo X. Miqueles - Elias S.Helou - Nikolay Koshev - Rafael F.C.Vescovi - Joao C. Cerqueira',
    author_email='eduardo.miqueles@lnls.br',
    
    description='Reconstructions Algorithms For Tomography',
    keywords=['tomography', 'reconstruction', 'imaging'],
    url='http://www.',
    download_url='',
    
    license='BSD',
    platforms='Any',
    install_requires = install_requires,
    
    classifiers=['Development Status :: 4 - Beta',
                 'License :: OSI Approved :: BSD License',
                 'Intended Audience :: Science/Research',
                 'Intended Audience :: Education',
                 'Intended Audience :: Developers',
                 'Natural Language :: English',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3.0',
                 'Programming Language :: C',
                 'Programming Language :: C++']
    
)


