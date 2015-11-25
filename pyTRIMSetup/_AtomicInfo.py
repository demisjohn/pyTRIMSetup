# -*- coding: utf-8 -*-
"""
AtomicInfo
Part of pyTRIMSetup by Demis D. John, Nov. 2015

Contains information about various atoms, as reported by TRIM.IN
All this info was originally found in the TRIM/SRIM program by James Zeigler (srim.org)

You can find info on other elements by using SRIM to create a target containing those elements,
then start running the implant simulation & abort it.  The program creates a TRIM.IN file
in it's directory, which contains all the info on each element in the target - copy that
data into the python lists here.

@author: Demis D. John
@email: demis@praevium.com

"""
####################################################
# Module setup etc.

#from __future__ import division  # Fix nonsense division in Python2.x (where 1/2 = 0 )- unneeded in new python versions
#import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
#import scipy as sp  # SciPy (signal and image processing library)

#import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
#import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
#from pylab import *              # Matplotlib's pylab interface, to enable typing commands just like MatLab.
#plt.ion()                            # Turned on Matplotlib's interactive mode

## For LaTeX usage (slows down plotting by maybe 5-10 sec per plot):
#from matplotlib import rc
#rc('text', usetex=True)        
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})

####################################################

# Add other modules here...

####################################################

print 'Loading Atomic info...'

_els = ['H', 'Al', 'Ga', 'As', 'In', 'Sb']
_nums = [1, 13, 31, 33, 49, 51]
_masses = [1.008, 26.982, 69.72, 74.922, 114.82, 121.75]
_surfbinding = [1.0, 3.36, 2.82, 1.26, 2.49, 2.72]
_displacement = [25 for x in _els]    # set to constant 25
_binding = [3 for x in _els]   # set to constant 3



