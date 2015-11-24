# -*- coding: utf-8 -*-
"""

TRIM Setup Example Using pyTRIMSetup module, by Demis D. John, Nov. 2015

A Python script to generate a TRIM.IN file for running TRIM.exe by James Zeigler (srim.org).
Nov. 2015, Demis D. John, Praevium Research Inc.

---------------------------------------------------------------

Python/Spyder script
Run Configuration should either:
    Execute in Current Interpreter
        OR
        
    Execute in New Interpreter 
        AND
    Interact with shell after completion (for plt.show() to work)


@author: Demis John
@email: demis@praevium.com


"""
####################################################
# General Module setup etc.

import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
import scipy as sp  # SciPy (signal and image processing library)
import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax

####################################################

# Add other modules here...
from pyTRIMSetup import *   # import the pyTRIMSetup module

####################################################

print 'Running...'




options = {}
options['Title'] = 'Testing the script' 
options['NumIons'] = 5000
options['AutoSaveNum'] = options['NumIons']
options['SimType'] = 2      # 1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons
options['RandomSeed'] = 0
options['Reminders'] = 0
options['DiskFiles'] = [0,0,0,0,0,0]    # booleans for: Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file
options['PlotType'] = 5     # 5: no plot
options['PlotExtents'] = (0,0)      # Xmin, Xmax = 0,0 for automatic

GaSb = Material(['Ga','Sb'], [0.5,0.5], 3.657)
AlAs = Material(['Al','Ga','As'], [0.3,0.3,0.4], 2.758)
AlSb = Material(['Al','Sb'], [0.5,0.5], 1.97)
target = Stack(  GaSb(1500) + AlAs(750) + GaSb(2000) + AlSb(2500) )     # top to bottom

target.implant(    Ion('H', 10, 7)    )     # Ion(ElementAbbrv, Energy_keV, Angle_degrees)

target.output('TestOutput.in', options=options, overwrite=True)





print 'done.'



