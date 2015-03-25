# -*- coding: utf-8 -*-
try:
    import pkg_resources 
    __version__ = pkg_resources.require("pyraft")[0].version
except:
    pass


from .raftypes import *
from .backprojection import *
from .sinogram import *
from .inversion import *
from .iterative import *
from .filters import *
from .phantom import *

from .iterative import *

try:
	import matplotlib.pyplot as plt
	from .misc import *
except:
	print("Matplotlib not found, not using plot functions")


try:
	import mayavi
	from .viswrapper import *
except:
	print("Mayavi not found, not using render functions")



