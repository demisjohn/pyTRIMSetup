# -*- coding: utf-8 -*-
"""

TRIM Setup Example Using pyTRIMSetup module, by Demis D. John, Nov. 2015

A Python script to generate a TRIM.IN file for running the ion-implantation simulator 
    TRIM.exe by James Zeigler (http://srim.org).
Syntax based on CAMFR by Peter Beinstman (http://camfr.sourceforge.net)
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
#import scipy as sp  # SciPy (signal and image processing library)
#import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
#import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax

####################################################

# Add other modules here...
from pyTRIMSetup import *   # import the pyTRIMSetup module

####################################################

print 'Running Script...'

'''
Get help on available arguments etc. by typing:
    >>> help( Material )
    >>> help( Ion )    
    >>> help( Stack )
    >>> import pyTRIMSetup as pt    # instead of importing `as *`
    >>> help( pt )
Or list all available functions/variables of an object:
    >>> AlAs = Material(['Al','As'],  [50,50],      3.752)
    >>> dir( AlAs )
    >>> help( AlAs )        # same as `help( Material )`
'''



'''
Setup target materials

A Material object is created as so:
    NewMat = Material(  [ListOfElements],  [ListOfMoleRatios],  Density)
    Defaults for optional args: CompoundCorrect=1, IsGas=False
'''
GaAs = Material(['Ga','As'], [0.5,0.5],     5.320)
AlGaAs = Material(['Al','Ga','As'], [98,02,100],   3.717)   # Mole-ratios will be normalized automatically
AlAs = Material(['Al','As'],  [50,50],      3.752)

# You can add custom elements that aren't yet included in AtomicInfo.py, as so:
X = Element('X', atnum=10, mass=9.01, surface_binding=10.2, binding=22, displacement=11)
GaX = Material( [ 'Ga', X  ], [0.5,0.5], 7.26 )     # note the lack of quotes'' around X - it is an Element object!

'''
Setup Target layer Stack

The Stack's Layers are created from top-to-bottom, by adding materials together, 
    while providing a thickness (Angstroms) for each layer.
    The underlying structures are python Lists, so you can group(), multiply*, and add+ Layers accordingly:
'''
# One period of AlAs & GaAs is multipled by two to repeat it twice:
repeatingpart = 2 * (  AlGaAs(150) + GaAs(110)  )   
# Insert this into the main Stack:  
target = Stack(  GaAs(1500) + repeatingpart + AlAs(2500, name='n-contact')  )     # top to bottom
# Added custom layer name to last AlAs layer.  Default name is just the compound, eg "AlAs"


'''
Setup ion to implant

The Ion is defined as
    Ion(  ElementAbbrv, Energy_keV, Angle_degrees )
    Angle is optional and defaults to 0-degrees.
'''
target.implant(    Ion('H', 10, 7)    )     # Ion(ElementAbbrv, Energy_keV, Angle_degrees)


'''
Setup a dictionary of simulator options:
    Eventually these will have default values, but for now you need to specify all of them.
'''
options = {}
options['Title'] = 'Testing the script with AlGaAs' 
options['NumIons'] = 5000       # Max number of ions to simulate - simulator stops at this many ions
options['AutoSaveNum'] = options['NumIons']     # Save the simulation at this interval
options['SimType'] = 2      # Cascades: 1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons.  See SRIM Menu for more info.
options['RandomSeed'] = 0
options['Reminders'] = 0
options['DiskFiles'] = [0,0,0,0,0,0]    # booleans for: Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file
options['PlotType'] = 5     # 5: no plot.  See the SRIM Plots menu.
options['PlotExtents'] = (0,0)      # Xmin, Xmax = 0,0 for automatic


'''
Generate the output TRIM.IN file with specified path/filename, passing the above options dictionary.
'''
target.output('ExampleOutput.in', options=options, overwrite=True)

'''
The generated file can now be copied to TRIM.IN in the SRIM/TRIM directory.
If you just run TRIM.exe, it will pick up all the simulation settings from this file.
'''


'''
Extraneous functionality:
'''
# Get more info on an object:
print target
# Get info on an Element:
print Element('Si')


print 'done.'



