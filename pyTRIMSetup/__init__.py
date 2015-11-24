# -*- coding: utf-8 -*-
"""

py TRIM Setup

A Python module for generating a TRIM.IN file for running TRIM.exe by James Zeigler (srim.org).
SRIM/TRIM is an ion-implantation monte-carlo simulator.

Nov. 2015, Demis D. John, Praevium Research Inc.

---------------------------------------------------------------
Example:

    import pyTRIMSetup as pt
    
    # Setup global options dictionary
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
    
    # Setup the ion to implant:
    Imp = pt.Ion('H', 10, 7)        # Ion(ElementName, Energy_keV, Angle_degrees)
    
    # Setup the target to bombard
    GaSb = pt.Material(['Ga','Sb'], [0.5,0.5], 3.657)   # ListOfElements, ListOfMoleFractions, Density
    AlAs = pt.Material(['Al','Ga','As'], [0.3,0.3,0.4], 2.758)
    AlSb = pt.Material(['Al','Sb'], [0.5,0.5], 1.97)
    target = pt.Stack(  GaSb(1500) + AlAs(750) + GaSb(2000) + AlSb(2500) )     # Tthicknesses, from top to bottom

    # Write the .IN file
    target.output('TestOutput.in', options=options, overwrite=True)



"""


print 'Loading pyTRIMSetup module...'

####################################################
# Module setup etc.
#from __future__ import division  # Fix nonsense division in Python2.x (where 1/2 = 0 )- unneeded in new python versions
import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
#import scipy as sp  # SciPy (signal and image processing library)
import os
import sys
from time import strftime

####################################################

# Internal module setup

ptDEBUG = False   # enable debugging output?

import _AtomicInfo as atom   # Conatins info on each atom, as found in SRIM

templatefile = 'TRIM.IN - Template'

templates_dir = os.path.join(os.path.dirname(_AtomicInfo.__file__), 'templates')
templatefile = os.path.join(templates_dir, templatefile)

'''
import pkg_resources
resource_package = atom.__name__  ## Could be any module/package name.
resource_path = os.path.join('/', templatefile)
print resource_package, resource_path
templatefile = pkg_resources.resource_string(resource_package, resource_path)
'''
####################################################

class Element(object):
    ''' Define an Element by passing the Abbreviation eg. 'H', 'Si' etc.
        
        Ga = Element('Ga')
        
        optional args:
            displacement (float): displacement energy
            binding (float): binding energy
            surface_binding (float): surface binding energy
    '''
    def __init__(self, *args, **kwargs):
        
        if len(args) >= 1:
            self.name, self.elnum, self.mass, self.surfbinding, self.binding, self.displacement = \
                self.element_lookup( str(args[0]) )
            self.displacement = kwargs.pop('displacement', self.displacement)
            self.binding = kwargs.pop('binding', self.binding)
            self.surfbinding = kwargs.pop('surface_binding', self.surfbinding)
        
        if len(args) >= 2:
            raise ValueError("Too many arguments!")
        
    #end __init__
    
    def element_lookup(self,  s, ion_only=False):
        '''Returns atomic info for an element, specified by it's abbreviations eg. 'H', 'Si', 'C' etc.
        returns ElementObject, AtomicNumber, AtomicMass, SurfaceBindingEnergy, LatticeBindingEnergy, DisplacementEnergy
        '''
        I = np.where(  np.array([s]) == np.array(atom._els)  )[0][0]
        
        return atom._els[I], atom._nums[I], atom._masses[I], atom._surfbinding[I], atom._binding[I], atom._displacement[I]
        
#end class(elements)



class Ion(Element):
    '''Class to define the ion being implanted.
        # implant Hydrogen at 150keV, 7° angle:
        H_implant = ion('H', 150, 7) 
        '''
    def __init__(self, *args):
        super(Ion, self).__init__()     # init `element` class, to get element lists
        if len(args)>=1:
            el = args[0]
            super(Ion, self).__init__(el)     # init `element` class, to get element lists
            #self.name, self.elnum, self.mass = self.element_lookup( str(el), ion_only=True )
            #if ptDEBUG: print self.element_lookup( str(el) )
        else:
            self.el, self.name, self.elnum, self.mass = None, None, None, None
        
        if len(args)>=2:
            self.energy = float( args[1] )
            if ptDEBUG: print self.energy
        else:
            self.energy = None
        
        if len(args)>=3:
            self.angle = float( args[2] )
            if ptDEBUG: print self.angle
        else:
            self.angle=0
        
        if len(args)>=4:
            raise ValueError("Too many arguments passed, max 3 arguments.")
    #end __init__
#end class(ion)

class Material(Element):
    '''Create a target material.
        Pass compounds as so:
            newmat = Material(  [list,of,elements], [list,of,fractions], Density, [CompoundCorr=1, Gas_Boolean=False])
            GaAs = Material( ['Ga', 'As'], [0.5, 0.5], 3.625)
            Al = Material( ['Al'], [1.0] )
        Optional arg name='p-contact' - if omitted, element names will be used
    '''
    def __init__(self, *args, **kwargs):
        super(Material, self).__init__()     # init `element` class, to get element lists
        if len(args)>=1:
            els = args[0]
            self.element, self.name, self.elnum, self.mass = [], [], [], []
            for el in els:
                elmt = Element(el)
                #name, elnum, mass = self.element_lookup( str(el) )
                if ptDEBUG: print self.element_lookup( str(el) )
                self.name.append(elmt.name)
                self.elnum.append(elmt.elnum)
                self.mass.append(elmt.mass)
                self.element.append(elmt)   # save the element object too
            self.description = kwargs.pop('name', None)
        else:
            self.el, self.name, self.elnum, self.mass = None, None, None, None
        
        if self.description == None:
            self.description = ''
            for n in self.name:
                self.description += n       # make description out of passed elements
        
        if len(args)>=2:
            if len(args[1]) != len(self.name): 
                    raise ValueError("Number of elements in 1st args & 2nd arg must match - need exactly one Mole Ratio for each Element provided.")
            
            self.molefrac = [float(x)/np.sum(args[1]) for x in args[1]] # check if can convert to number & normalize
            if ptDEBUG: print self.molefrac
        else:
            self.molefrac = None
        
        '''
        if len(args)>=3:
            self.thickness = float( args[2] )
        else:
            self.thickness = None
        '''
        
        if len(args)>=3:
            self.density = float( args[2] )
        else:
            self.density = None     # replace with automatic interpolated calculation?
        
        if len(args)>=4:
            self.compoundCorrection = float( args[3] )
        else:
            self.compoundCorrection = 1.0
        
        if len(args)>=5:
            self.isGas = bool( args[4] )
        else:
            self.isGas = False
    
    def __add__(self,other):
        return [self, other]
    
    def __call__(self, thickness):
        return [Layer(self, thickness)]     # return list, so can add Layers later
    
#end class material


class Layer(object):
    '''Invisible to user.  Add thickness to a material.
        Layer(Material_Object, thickness_angstroms)
        '''
    
    def __init__(self, MaterialObj, thickness):
        if not isinstance(MaterialObj, Material):
            raise ValueError("First argument should be a Material object!")
        
        # copy Material's attributes
        self.material = MaterialObj
        self.name = MaterialObj.name
        self.description = MaterialObj.description
        self.elnum = MaterialObj.elnum
        self.mass = MaterialObj.mass
        self.molefrac = MaterialObj.molefrac
        self.thickness = thickness
        self.density = MaterialObj.density
        self.compoundCorrection = MaterialObj.compoundCorrection
        self.isGas = MaterialObj.isGas

    def __add__(self,other):
        '''addition: concatenate to list'''
        return [self, other]
        
        

class Stack(object):
    ''' stack multiple target materials up, from top to bottom.
        Pass a list, like so:
            Target = Stack(  GaAs + InAs + InAs_thick )  '''
    def __init__(self, *args):
        if ptDEBUG: print "Stack.init():  ", args[0]
        self.stack = args[0]
        
        ## Generate elements list
        elnames, els = [], []
        for s in self.stack:
            # each layer in the Stack
            for e in range(   len(s.material.name)   ):
                # each element in the Layer
                elnames.append( s.material.name[e] )
                els.append( Element(elnames[-1]) )  # make new element object from this one's name
        #elnames_u, elnames_i = np.unique(elnames, return_index=True)    # will sort the elements
        #if ptDEBUG: print els
        #self.elements = [els[i]   for i in elnames_i]  # save unique element objects
        self.elements = els  # save all elements, even if repeated
        
        for i, e in enumerate( self.elements ):
            e.idnum = i+1       # ID number of this element
            
    #end __init__
    
    
    def __len__(self):
        return len(self.stack)
    
    def __repr__(self):
        return str(self.stack)
    
    def __add__(self,other):
        return [self, other]
    
    def get_numElements(self):
        '''Return number of unique elements contained.'''
        return len(self.elements)
    
    def implant(self, ion_object):
        '''Define the Ion to implant.  Takes a single Ion object as input.
        Use as so:
            Target.implant(  Ion('H', 10, 7)  )
        '''
        self.ion = ion_object
    
    def output(self, filepath, options=None, overwrite=False, warn=True):
        '''Write the *.IN output file to `filepath`.  The .IN extension will be added automatically.
            To use the file in TRIM.exe, the file should be renamed to TRIM.IN and placed in the SRIM/TRIM folder.
            
            options : dictionary containing global options.
                Example:
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
                    
            
            overwrite : {True | False}
                Overwrite existing files? False by default.
            
            warn : {True | False}
                Issue warning when overwriting a file?  True by default.
        '''
        if os.path.exists(filepath):        
            if not overwrite:
                raise IOError( "File `%s` already exists, aborting.  Set `overwrite=True` to overwrite the file."%(filepath)  )
            else:
                if warn:
                    print "WARNING: Overwriting file at: \n\t\%s"%(filepath)
        #end if(file-exists)
        

        le = '\r\n'   # line-ending: put this at the end of every line!
        tab = '    '    
        
        ################################
        # Start writing the output file
        ################################
        
        # open TRIM.IN template file:
        if not os.path.exists(templatefile):        
            raise IOError( "Could not find the TRIM.IN template file in the module directory!  Looked at path:\n\t\%s" % (templatefile)  )
        t = open(templatefile, 'r')
        ts = t.readlines();     t.close()
        
        
        
        
        f = open(filepath,'w')
        f.write('==> SRIM-2013.00.  Generated by pyTRIM, Demis D. John 2015.' + le) # 1st comment line
        f.write(ts[1])  # comment line: Ion: Zi...
        f.write(tab + str(self.ion.elnum) + tab + str(self.ion.mass) + tab + str(self.ion.energy) + tab + \
             str(self.ion.angle) + tab + str(options['NumIons']) + tab + str(1) + tab + str(options['AutoSaveNum']) + le     )
        
        f.write(ts[3])  # comment line: Cascades...
        f.write(tab*3 + str(options['SimType']) + tab*3 + str(options['RandomSeed']) + tab + str(options['Reminders']) + le)
        
        f.write(ts[5])  # comment: Diskfiles...
        fstr = tab*2
        for d in options['DiskFiles']:
            fstr += str(d) + tab
        f.write(fstr + le)
        
        f.write(ts[7])  # comment: Target Material...
        f.write('"%s"' %(options['Title'])  +  tab + str(self.get_numElements()) + tab + str(len(self)) + le)
        
        f.write(ts[9])  # comment: PlotType...
        f.write(tab + str(options['PlotType']) + tab + str(options['PlotExtents'][0]) + tab + str(options['PlotExtents'][1]) + tab + le)
        
        ## List of Elements:
        f.write(ts[11])  # comment: Target Elements...
        for e in self.elements:
            f.write("Atom %i = %s =      %f  %f" %(e.idnum, e.name, e.elnum, e.mass)    + le)
        
        ## List of Target layers:
        f.write(ts[13])  # comment: Layer...
        f.write(ts[14])  # comment: Numb.  ...
        # 1      "GaSb"           4830  6.294      .5      .5       0       0       0
        prevelt = 0 # position of previous element
        for i,l in enumerate(self.stack):
            numelts = 0  # current position
            
            fstr = ' %i    "%s"'%(i+1, l.description) + tab + str(l.thickness) + tab + str(l.density) + tab
            
            for ii in range(prevelt):
                # write 0's for mole frac's on elements not in this layer
                numelts = numelts + 1
                fstr += '0' + tab
                
            for n, m in enumerate(l.molefrac):
                numelts = numelts + 1
                prevelt = numelts       # record position of last molefrac entered
                fstr += str(m) + tab
            
            for ii in range(self.get_numElements() - numelts ):
                # fill rest of elements with 0 mole frac.
                fstr += '0' + tab
                
            fstr += le    
            f.write(fstr)
        #end for(Layers)
        
        f.write(ts[16])  # comment: Target layer phases...
        f.write(' ')
        for l in self.stack:
            B = 1   if  l.isGas else  0
            f.write( str(B) + ' ' )
        f.write( le )
        
        f.write(ts[18])  # comment: Target compound corrections...
        f.write(' ')
        for l in self.stack:
            f.write( str(l.compoundCorrection) + '   ' )
        f.write( le )
        
        f.write(ts[20])  # comment: atom displacements...
        for e in self.elements:
            f.write( tab + str(e.displacement)  )
        f.write( le )
        
        f.write(ts[22])  # comment: atom lattice binding...
        for e in self.elements:
            f.write( tab + str(e.binding)  )
        f.write( le )
        
        f.write(ts[24])  # comment: atom surface binding...
        for e in self.elements:
            f.write( tab + str(e.surfbinding)  )
        f.write( le )
        
        f.write(ts[26])  # comment: Stopping power verions
        f.write(ts[27])  # Stopping power verions
        
        f.close()
            
#end class(Stack)
    
        


