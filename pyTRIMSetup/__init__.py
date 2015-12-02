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
from time import strftime   # get date/time as string

####################################################

# Internal module setup

global ptDEBUG
ptDEBUG = False   # enable debugging output?

import _AtomicInfo as atom   # Conatins info on each atom, as found in SRIM

_SRIMMaxLayers = 51     # Maximum no. layers SRIM/TRIM can take


from sys import platform as _platform       # Used in Stack.output



''' This is now hard-coded
templatefile = 'TRIM.IN - Template'
templates_dir = os.path.join(os.path.dirname(_AtomicInfo.__file__), 'templates')
templatefile = os.path.join(templates_dir, templatefile)
'''

'''
import pkg_resources
resource_package = atom.__name__  ## Could be any module/package name.
resource_path = os.path.join('/', templatefile)
print resource_package, resource_path
templatefile = pkg_resources.resource_string(resource_package, resource_path)
'''
####################################################

#  These override the value set above in `ptDEBUG`
def set_DEBUG():
    '''Enable verbose output for debugging.'''
    global ptDEBUG
    ptDEBUG = True

def unset_DEBUG():
    '''Disable verbose debugging output.'''
    global ptDEBUG
    ptDEBUG = False

def DEBUG():
    '''Returns whether DEBUG is true or false'''
    return ptDEBUG


def importElements( trimfile ):
    '''
    Use this function to import elements from a TRIM.IN file created by SRIM.exe.
    The Elements stored in the file will be returned in a list, containing the corresponding binding energies etc.
    
    Parameters
    ----------
    
    trimfile : path
        Path to TRIM.IN file you'd like to import elements from.
    
    
    Returns
    -------
    
    Elements : list of Element objects
        Lists of Element objects, which contains the surface binding energies, lattice binding energies and displacement energies, element numbers and atomic masses.
        This list can be passed directly to `exportElements()` to create a new _AtomicInfo.py file.
    
    '''
    
    import re   # RegEx
    
    if not os.path.exists(trimfile):        
        raise IOError( "Could not find the specified TRIM.IN!  Looked at path:\n\t\%s" % (trimfile)  )
    t = open(trimfile, 'r')
    ts = t.readlines();     t.close()
    
    '''
    # Find first word:
    atompat = re.compile(     r'\s*(\w*)\s*.*' ,    flags=(re.IGNORECASE) )
    '''
    
    # Find Atom data
    pat = re.compile(     r'\s*Atom\s*(?P<pos>\d*)\s*=\s*(?P<name>\w*)\s*=\s*(?P<elnum>\d+\.?\d*)\s*(?P<mass>\d+\.?\d*)' ,    flags=(re.IGNORECASE) )
    
    
    cont=True
    i=11    # start reading at 1st line of target elements
    pos, name, elnum, mass = [],[],[],[]
    while cont:
        i+=1
        L = ts[i]   # current line of file
        
        # perform the search:
        m = pat.search(  L  )      # search for RegEx pattern
        # m will contain any 'groups' () defined in the RegEx pattern.
        if m:
            # found regex pattern on this line
            #(?P<pos>\d*)  (?P<name>\w*)  (?P<elnum>\d+\.?\d*)  (?P<mass>\d+\.?\d*)
            pos.append(   int(  m.group('pos')  )    )	# grab named group from RegEx
            name.append(    m.group('name')    )
            elnum.append(    float(  m.group('elnum')  )    )
            mass.append(    float(  m.group('mass')  )    )
            if DEBUG(): print 'Atom found on line#%i:'%(i), m.groups(), "\n\tLast (pos,name,elnum,mass) = ",pos[-1],name[-1],elnum[-1],mass[-1]
            
        #groups() prints all captured groups
        else:
            # did not find the regex pattern on this line
            if DEBUG(): print "No Atom found at line #%i.  "%( i )
            lastline =  i
            cont=False  # unnecessary if `break` works
        #raw_input()
    #end while(atoms)
    
    # Find Layer rows:
    layernum = re.compile(     r'\s*(\d+)\s*"(.*)"\s*(\d+\.?\d*)\s*(\d+\.?\d*).*' ,    flags=(re.IGNORECASE) )
    #  10     "Layer 10"           10000  17.88402     0.5000
    
    i = lastline+1
    cont=True
    while cont:
        i+=1
        L = ts[i]   # current line of file
        #if DEBUG(): print "Layer search: Line %i:\n\t"%(i), L
        
        # perform the search:
        m = layernum.search(  L  )      # search for RegEx pattern
        # m will contain any 'groups' () defined in the RegEx pattern.
        if m:
            d = m.group(1)
            if DEBUG(): print 'Layer found on line#%i:'%(i), m.groups()
        #groups() prints all captured groups
        else:
            if DEBUG(): print "No Layer found at line #%i."%( i )
            lastline =  i
            cont=False  # unnecessary if `break` works
            break   # end the while loop
    #end while(target layers)
    #if DEBUG(): print "lastline=",lastline
    
    # Get displacement E's:
    disp = ts[lastline + 5]
    #if DEBUG(): print "disp=", disp
    disp = disp.split()     # convert to list, splitting at spaces
    disp = [float(x) for x in disp]     # convert str --> float
    
    # Lattice binding E's
    latbind = ts[lastline + 7]
    #if DEBUG(): print "latbind=", latbind
    latbind = latbind.split()
    latbind = [float(x) for x in latbind]
    
    # Surface binding E's
    surfbind = ts[lastline + 9]
    #if DEBUG(): print "surfbind=", surfbind
    surfbind = surfbind.split()
    surfbind = [float(x) for x in surfbind]
    
    #old: Elements = [pos,name,elnum,mass]
    Elements=[]
    for ii in range(  len(name)  ):
        # Generate new Element objects for each of these & add to list:
        Elements.append(  Element(name=name[ii], atnum=elnum[ii], mass=mass[ii], \
            surface_binding=surfbind[ii], lattice_binding=latbind[ii], displacement=disp[ii])  )
    #end for(name)
    
    #if DEBUG(): print [str(x) for x in Elements]
    return Elements     #old:   Elements, surfbind, latbind, disp
#end importElements()        


def exportElements( elements, filepath, overwrite=False, warn=True ):
    '''Export a list of Element objects to an AtomicInfo.py file, for performing element lookups with.  The resulting file is a Python script, which pyTRIMSetup can import as part of the module, typically stored in pyTRIMSetup/_AtomicInfo.py.
    
    Parameters
    ----------
    elements : list
        A list of Element objects.  The output from `importElements()` can be passed directy to this parameter.
    
    filepath : path (str)
        Where to save the output file to.  The output file should be saved to pyTRIMSetup/_AtomicInfo.py to allow Element lookups form the new file.
    
    overwrite : {True | False}, optional
        Overwrite existing files? False by default.
    
    warn : {True | False}, optional
        Issue warning when overwriting a file?  True by default.
    '''
    
    # Check for file existence:
    if os.path.exists(filepath):        
        if not overwrite:
            raise IOError( "File `%s` already exists, aborting.  Set `overwrite=True` to overwrite the file."%(filepath)  )
        else:
            if warn:
                print "WARNING: Overwriting file at: \n\t\%s"%(filepath)
    #end if(file-exists)

    print "Writing to file: %s" %(filepath)
    f = open(filepath,'w')
    fstr = '# -*- coding: utf-8 -*-' + '\n'
    fstr += '"""\nAtomicInfo\n' + '\n'
    fstr += '\tGenerated on: %s\n' % (  strftime("%Y-%m-%d %H.%M.%S")  )
    fstr += '\tCreated by pyTRIMSetup.exportElements(), Demis D. John, 2015\n'
    fstr += '\tContains information about various atoms, as reported by TRIM.IN\n'
    fstr += '\tAll this info was originally found in the TRIM/SRIM program by James F. Zeigler (srim.org).\n'
    fstr += '"""\n\n'
    fstr += 'print "Loading Atomic info..."\n\n'
    f.write(fstr) 
    
    # element names:
    fstr = "_els = %s ; \n\n" %(   str( [x.name  for x in elements] )   )
    f.write( fstr )
    
    # element number
    fstr = "_nums = %s ; \n\n" %(   str( [x.elnum  for x in elements] )   )
    f.write( fstr )
    
    # element mass
    fstr = "_masses = %s ; \n\n" %(   str( [x.mass  for x in elements] )   )
    f.write( fstr )
    
    # displacement E
    fstr = "_displacement = %s ; \n\n" %(   str( [x.displacement  for x in elements] )   )
    f.write( fstr )
    
    # lat binding E
    fstr = "_binding = %s ; \n\n" %(   str( [x.binding  for x in elements] )   )
    f.write( fstr )
    
    # surf binding E
    fstr = "_surfbinding = %s ; \n\n" %(   str( [x.surfbinding  for x in elements] )   )
    f.write( fstr )
    
    f.close()
#end exportElements() 


class Element(object):
    ''' Define an Element by passing the Abbreviation eg. 'H', 'Si' etc.
        
        >>> Silicon = Element('Si')      # will perform a lookup of the Element to obtain parameters
        
        Keyworded arguments can be used to override the default parameters for an element, eg.
        
        >>> Silicon = Element('Si', mass=28.05)     # Only mass is overidden
        
        If only keyword-args are passed, then no lookup will be performed - it is assumed you are defining the Element manually (and all parameters will thus need to be defined by you).
        
        Manual Element definition (provide a `name=` parameter):
        
        >>> demisium = Element( name='De', atnum=149, mass=297.997, displacement=25, lattice_binding=3, surface_binding=1.097)
        
        
        optional args:
            name (str) : name for the element
            atnum (int) : atomic number
            mass (float) : atomic mass in a.m.u.
            displacement (float): displacement energy (keV)
            lattice_binding (float): binding energy (keV)
            surface_binding (float): surface binding energy (keV)
    '''
    def __init__(self, ElementName=None, **kwargs):
        
        if ElementName:
            if DEBUG(): print "Got ElementName: %s"%ElementName
            name, self.elnum, self.mass, \
            self.surfbinding, self.binding, self.displacement = \
                self.element_lookup( str(ElementName) )
            
            if name==None:
                self.name = ElementName
            else:
                self.name = name
            
            # Optional keyword args
            self.elnum = kwargs.pop('atnum', self.elnum)
            self.mass =  kwargs.pop('mass', self.mass)
            self.displacement = kwargs.pop('displacement', self.displacement)
            self.binding = kwargs.pop('lattice_binding', self.binding)
            self.surfbinding = kwargs.pop('surface_binding', self.surfbinding)
            
        else:
            # no ElementName passed, use all kwargs to define.
            # Optional keyword args
            self.name = kwargs.pop('name', None)
            self.elnum = kwargs.pop('atnum', None)
            self.mass =  kwargs.pop('mass', None)
            self.displacement = kwargs.pop('displacement', None)
            self.binding = kwargs.pop('lattice_binding', None)
            self.surfbinding = kwargs.pop('surface_binding', None)
        #end if(ElementName)
            
        # Make sure element is defined!
        if self.name==None or self.elnum==None or self.mass==None or self.surfbinding==None or self.binding==None or self.displacement==None:
            ErrStr = "Please provide all parameters for your custom Element! Found these parameters:\n\tname= " + str(ElementName) + "\n\tatnum= " + str(self.elnum) + "\n\tmass= " + str(self.mass) + "\n\tsurface_binding= " + str(self.surfbinding) + "\n\tbinding= " + str(self.binding) + "\n\tdisplacement= " + str(self.displacement)
            raise ValueError(ErrStr)
        
    #end __init__
    
    def __str__(self):
        '''How to `print` this object.'''
        '''
        Abbreviation: 'H'
        Atomic Number: 1
        Atomic Mass: 1.008
        Surface Binding E: 30 keV
        etc. etc.
        '''
        pstr = ""   #"<pyTRIMSetup Element Object>\n"
        pstr += "Abbreviation: '%s'\n" % (self.name)
        pstr += "Atomic Number: %i\n" % (self.elnum)
        pstr += "Atomic Mass: %f amu\n" % (self.mass)
        pstr += "Surface Binding Energy: %f keV\n" % (self.surfbinding)
        pstr += "Lattice Binding Energy: %f keV\n" % (self.binding)
        pstr += "Displacement Energy: %f keV" % (self.displacement)
        return pstr
    
    def element_lookup(self,  s, ion_only=False):
        '''Returns atomic info for an element, specified by it's abbreviations eg. 'H', 'Si', 'C' etc.
        returns ElementObject, AtomicNumber, AtomicMass, SurfaceBindingEnergy, LatticeBindingEnergy, DisplacementEnergy
        '''
        try:
            I = np.where(  np.array([s]) == np.array(atom._els)  )[0][0]
        except IndexError:
            return [None]*6
        return atom._els[I], atom._nums[I], atom._masses[I], \
            atom._surfbinding[I], atom._binding[I], atom._displacement[I]
        
#end class(elements)



class Ion(Element):
    '''Class to define the ion being implanted.
        # implant Hydrogen at 150keV, 7° angle:
        H_implant = ion('H', 150, 7) 
        '''
    def __init__(self, *args):
        #super(Ion, self).__init__()     # init `element` class, to get element lists
        if len(args)>=1:
            el = args[0]
            super(Ion, self).__init__(el)     # init `element` class, to get element lists
            
            #self.name, self.elnum, self.mass = self.element_lookup( str(el), ion_only=True )
            #if DEBUG(): print self.element_lookup( str(el) )
        else:
            self.el, self.name, self.elnum, self.mass = None, None, None, None
        
        if len(args)>=2:
            self.energy = float( args[1] )
            if DEBUG(): print self.energy
        else:
            self.energy = None
        
        if len(args)>=3:
            self.angle = float( args[2] )
            if DEBUG(): print self.angle
        else:
            self.angle=0
        
        if len(args)>=4:
            raise ValueError("Too many arguments passed, max 3 arguments.")
    #end __init__
    
    def __str__(self):
        '''How to `print` this object.'''
        '''
        Abbreviation: 'H'
        Atomic Number: 1
        Atomic Mass: 1.008
        etc. etc.
        '''
        pstr = ""   #"<pyTRIMSetup Ion Object>\n"
        pstr += "Abbreviation: '%s'\n" % (self.name)
        pstr += "Atomic Number: %i\n" % (self.elnum)
        pstr += "Atomic Mass: %f amu\n" % (self.mass)
        pstr += "Energy: %f keV\n" %(self.energy)
        pstr += u"Incidence Angle: %f°" %(self.angle)
        return pstr
    #end __str__
#end class(ion)

class Material(Element):
    '''Create a target material.
    
        Pass compounds as so:
        
            newmat = Material(  [list,of,elements], [list,of,fractions], Density, [CompoundCorr=1, Gas_Boolean=False])
            
        Example:
            GaAs = Material( ['Ga', 'As'], [0.5, 0.5], 3.456)
            Al = Material( ['Al'], [1.0], 4.567 )
        
        To make your own elements (that aren't present in the AtomicInfo.py lookup table), make an Elelemnt first and pass that:
            
        Example:
            X = Element('X', atnum=10, mass=9.01, surface_binding=10.2, binding=22, displacement=11)
            GaX = Material( [ 'Ga', X  ], [0.5,0.5], 7.26 )
            # Note the missing ''quotes around Element X - it is an Element object
        
        Optional argument:
            name : str
                Optional name for this material, eg. 'p-contact'.  This will become the layer name in TRIM.  If omitted, the constituent element names will be used.  
        
        Element abbreviations (eg. 'H', 'Si' etc.) are matched up to their parameters via an `element_lookup()` function, which looks in pyTRIMSetup/_AtomicInfo.py to find the info on that element. The functions importElements()/exportElements() can be used to create the _AtomicInfo.py file.
    '''
    def __init__(self, *args, **kwargs):
        #super(Material, self).__init__()     # init `element` class, to get element lists
        
        # 1st arg: element abbreviation
        if len(args)>=1:
            if isinstance(args[0], (tuple,list) ):
                els = args[0]
            elif isinstance(args[0], str):
                els = [  args[0]  ]    # add to list, so can iterate
            else:
                ErrStr = "Expected 1st argument to be list/tuple or string.  Instead got type: %s = `%s`."%( type(args[0]),  args[0] )
                raise ValueError( ErrStr )
            #end if(list or str)
            
            self.element, self.name, self.elnum, self.mass = [], [], [], []
            for el in els:
                if isinstance(el, Element):
                    elmt = el       # it is already an element object
                elif isinstance(el, str):
                    elmt = Element(el)  # convert to Element by passing name
                else:
                    ErrStr = "Unknown type passed as element:\n\t%s"%(type(el)) + "\nfor argument:\n\t%s"%(el)
                    raise ValueError(ErrStr)
                
                #name, elnum, mass = self.element_lookup( str(el) )
                if DEBUG(): print self.element_lookup( str(el) )
                self.name.append(elmt.name)
                self.elnum.append(elmt.elnum)
                self.mass.append(elmt.mass)
                self.element.append(elmt)   # save the element object too
            self.description = kwargs.pop('name', None)
        else:
            self.el, self.name, self.elnum, self.mass = None, None, None, None
        
        if self.description == None:
            self.description = ''
            for e in self.element:
                self.description += e.name       # make description out of passed element's names
        
        # 2nd arg: mole-frac.
        if len(args)>=2:
            if isinstance(args[1],  (tuple, list) ):
                mf = args[1]
            else:
                mf = [  args[1]  ]    # add to list, so can iterate
            #end if(list or str)
            
            if len(mf) != len(self.name): 
                    raise ValueError("Number of elements in 1st args & 2nd arg must match - need exactly one Mole Ratio for each Element provided.")
            
            self.molefrac = [   float(x)/np.sum(mf)   for   x in mf   ] # check if can convert to number & normalize
            if DEBUG(): print self.molefrac
        else:
            self.molefrac = None
        
        # 3rd arg: density
        if len(args)>=3:
            self.density = float( args[2] )
        else:
            self.density = None     # replace with automatic interpolated calculation?
        
        if len(args)<3:
            raise ValueError("Material() requires at least three arguments: ListOfElements, ListOfMoleRatios, Density")
        
        if len(args)>=4:
            self.compoundCorrection = float( args[3] )
        else:
            self.compoundCorrection = 1.0
        
        if len(args)>=5:
            self.isGas = bool( args[4] )
        else:
            self.isGas = False
        
        # get kwargs:
        self.compoundCorrection = kwargs.pop('CompoundCorrection', self.compoundCorrection)
        self.isGas = kwargs.pop('IsGas', self.isGas)
    #end init()
    
    def __str__(self):
        '''How to `print` this object.'''
        '''
        Al(0.5)As(0.5)
        '''
        pstr = ""   #"<pyTRIMSetup Material Object>\n"
        for e in range(len(self.element)):
            pstr+= "%s(%f)" % (self.element[e].name, self.molefrac[e])
        return pstr
    #end __str__
    
    def __add__(self,other):
        return [self, other]
    
    def __call__(self, thickness, name=None):
        '''Return a list containing Layer object with given thickness.  List allows __add__'ing later.  Optional `name=` param will be passed to Layer, overiding Material name.'''
        return [Layer(self, thickness, name=name)]     
    
#end class material


class Layer(object):
    '''Invisible to user.  Add thickness to a material.
        Layer(Material_Object, thickness_angstroms)
        
        The Material object returns a Layer object when called() with a thickness value.  The resulting Layer is the same as a Material but with an added thickness attribute.
        
        Optional `name=` param can be passed. If omitted, name is taken from the Material object.
        '''
    
    def __init__(self, MaterialObj, thickness, name=None):
        if not isinstance(MaterialObj, Material):
            raise ValueError("First argument should be a Material object!")
        
        # copy Material's attributes
        self.material = MaterialObj
        if name != None:
            self.description = name    # use passed name
        else:
            self.description = MaterialObj.description    # use Material's name
        self.elnum = MaterialObj.elnum
        self.mass = MaterialObj.mass
        self.molefrac = MaterialObj.molefrac
        self.thickness = thickness
        self.density = MaterialObj.density
        self.compoundCorrection = MaterialObj.compoundCorrection
        self.isGas = MaterialObj.isGas
    
    def __str__(self):
        '''How to `print` this object.'''
        '''
        [Name]: Ga(0.5)As(0.5) = 1500 Ang.
        '''
        pstr = ""   #"<pyTRIMSetup Layer Object>\n"
        pstr += u"'%s': %s = %f Ang." % (self.description, self.material, self.thickness)
        return pstr
    #end __str__
    
    def __add__(self,other):
        '''addition: concatenate to list'''
        return [self, other]
        
        

class Stack(object):
    ''' Stack multiple target materials up, from top to bottom.
    
        Pass Materials with some thickness (Angstroms), like so:
            Target = Stack(  GaAs_top(150) + InAs(250) + InAs_bottom(10000) )  
        
        Internally, the Material object returns a Layer object when called() with a thickness value.  Layer objects return a concatenated [List] of Layers when added together, and the Stack consists of this List of Layers.'''
    def __init__(self, *args):
        if DEBUG(): print "Stack.init():  ", args[0]
        self.stack = args[0]
        
        ## Generate elements list
        elnames, els = [], []
        for s in self.stack:
            # for each layer in the Stack
            for e in range(   len(s.material.element)   ):
                # for each element in the Layer
                elnames.append( s.material.element[e].name )
                els.append( s.material.element[e] )
                
                
                '''
                elnames.append( s.material.name[e] )
                els.append( Element(elnames[-1]) )  # make new element object from this one's name
                '''
        #end for(StackLayers/Materials)
        
        """
        ## Find only unique elements - don't repeat identical elements
        #   This should be re-added as an option later.  Currently don't have the login in place to order the elements properly
        elnames_u, elnames_i = np.unique(elnames, return_index=True)    # will sort the elements
        if DEBUG(): print els
        self.elements = [els[i]   for i in elnames_i]  # save unique element objects
        """
        self.elements = els  # save all elements, even if repeated
        
        self.elidnum = range(  1, len(self.elements)+1  )
        
    #end __init__
    
    
    def __len__(self):
        return len(self.stack)
    
    #def __repr__(self):
    #    return str(self.stack)
    # Let python handle __repr__, so you know what type of object it is etc.
    
    def __str__(self):
        '''How to `print` this object.'''
        '''
        Layer 1: Si(1.0)[1400Å]
        Layer 2: Ga(0.5)As(0.5)[1000Å]
        Total Thickness = 2400Å
        Number of Elements: 3
        '''
        pstr = ""   #"<pyTRIMSetup Stack Object>\n"
        for i,L in enumerate(self.stack):
            pstr += "Layer %i: %s\n" %(i+1, L)
        pstr += u"Total Thickness = %f Ang.\n" % ( self.get_thickness() )
        pstr += "Number of Elements: %i" %( self.get_numElements() )
        #if DEBUG(): print pstr
        return pstr
    #end __str__
    
    def __add__(self,other):
        '''Addition - return the list of layers to do math.'''
        return self.stack + list(other)
    
    def __radd__(self,other):
        '''Reverse Add - useful for adding to integers.
         - return the list of layers to do math.'''
        return list(other) + self.stack
    
    def __mul__(self, other):
        '''Multiply - return the list of layers to do math.'''
        return self.stack * int(other)
    
    def __rmul__(self, other):
        '''Reverse Multiply (useful for integer math) - return the list of layers to do math.'''
        return int(other) * self.stack
    
    def get_numElements(self):
        '''Return number of unique elements contained.'''
        return len(self.elements)
    
    def get_thickness(self):
        '''Return total thickness of this Stack (by adding up thicknesses of each contained layer).'''
        return np.sum(  [layer.thickness   for layer in   self.stack]  )
    
    def implant(self, ion_object):
        '''Define the Ion to implant.  Takes a single Ion object as input.
        Use as so:
            Target.implant(  Ion('H', 10, 7)  )
        '''
        self.ion = ion_object
    
    def output(self, filepath, options=None, overwrite=False, warn=True, split=False, numlayers=30):
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
                    
            
            overwrite : {True | False}, optional
                Overwrite existing files? False by default.
            
            warn : {True | False}, optional
                Issue warning when overwriting a file?  True by default.
            
            split : {True | False}, optional
                Split the output file up if the Stack contains more layers than SRIM/TRIM can handle?  (As of 2015, SRIM/TRIM gives an error if more than 51 layers are set.)  
                This allows you to simulate one implant, and then take the transmitted ions (TRANSMIT.txt - see SRIM manual Ch. 9-2) and launch those into the next section of the target.
                Output files will be numbered with 01 being at the top (first to be implanted).
            
            numlayers : int, optional
                If `split=True`: How many layers should the resulting files contain?  Defaults to 30.
        '''
        
        if not split:
            if len(self.stack) > _SRIMMaxLayers:
                print 'WARNING: Stack.output(): SRIM/TRIM.exe can not work with more than %i layers, and your structure contains %i layers!\nSRIM.exe/TRIM.exe may issue a "Runtime Error 9: Subscript out of range".  Use the `split=True` option to output multiple target files.'%( _SRIMMaxLayers, len(self.stack) )
        else:
            #SPLIT = True
            # set to output Transmitted ions:
            options['DiskFiles'][2] = 1     # Transmit.txt output
            #options['DiskFiles'][4] = 2     # Collisions.txt output
            
            # Reference: options['DiskFiles'] = [0,0,0,0,0,0]    # booleans for: Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file
        #end if(split)

        
        if split:
            allstacks = self.splitstack(numlayers=numlayers)  # return multiple stacks, split up according to _SRIMMaxLayers
        else:
            allstacks = [self]      # make list so it's iterable
        
        
        # This is now set at the top of the module.
        #le = '\r\n'   # (If run on MacOS) line-ending: put this at the end of every line!
        #le = '\r'       # If executed on WinXP
        tab = '    ' 

        # Determine line-ending based on system we're executing on. 
        #   _platform imported at top, from sys.platform  
        if _platform == "linux" or _platform == "linux2":
           # linux
           le = '\r\n'      # not tested!
        elif _platform == "darwin":
           # MAC OS X
           le = '\r\n'
        elif _platform == "win32":
           # Windows
           le = '\n'
        
        '''
        # open TRIM.IN template file:
        if not os.path.exists(templatefile):        
            raise IOError( "Could not find the TRIM.IN template file in the module directory!  Looked at path:\n\t\%s" % (templatefile)  )
        t = open(templatefile, 'r')
        ts = t.readlines();     t.close()
        '''
        
        tempfilepath = filepath
        
        for Is,stacks in enumerate(allstacks):
            '''iterate through each sub-stack'''
            
            ################################
            # Start writing the output file
            ################################
            
            if split:
                # construct new filename:
                if os.path.splitext(tempfilepath)[1].lower() == '.in'.lower():
                    # convert 'TRIM.IN' into 'TRIM 01.IN'
                    filepath = os.path.splitext(tempfilepath)[0] + "_%0.3i"%(Is) + os.path.splitext(tempfilepath)[1]
                else:
                    filepath = tempfilepath + "_%0.3i"%(Is)
            #end if(SPLIT)
            
            # Check for file existence:
            if os.path.exists(filepath):        
                if not overwrite:
                    raise IOError( "File `%s` already exists, aborting.  Set `overwrite=True` to overwrite the file."%(filepath)  )
                else:
                    if warn:
                        print "WARNING: Overwriting file at: \n\t\%s"%(filepath)
            #end if(file-exists)
        
            print "Writing to file: %s" %(filepath)
            f = open(filepath,'w')
            f.write('==> SRIM-2013.00.  Generated by pyTRIMSetup by Demis D. John 2015.' + le) # 1st comment line
            #f.write(ts[1])  # comment line: Ion: Zi...
            f.write("Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number."+le)
            f.write(tab + str(self.ion.elnum) + tab + str(self.ion.mass) + tab + str(self.ion.energy) + tab + \
                 str(self.ion.angle) + tab + str(options['NumIons']) + tab + str(1) + tab + str(options['AutoSaveNum']) + le     )
        
            #f.write(ts[3])  # comment line: Cascades...
            f.write("Cascades(1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons), Random Number Seed, Reminders"+le)
            f.write(tab*3 + str(options['SimType']) + tab*3 + str(options['RandomSeed']) + tab + str(options['Reminders']) + le)
        
            #f.write(ts[5])  # comment: Diskfiles...
            f.write("Diskfiles (0=no,1=yes): Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file"+le)
            fstr = tab*2
            for d in options['DiskFiles']:
                fstr += str(d) + tab
            f.write(fstr + le)
        
            #f.write(ts[7])  # comment: Target Material...
            f.write("Target material : Number of Elements & Layers"+le)
            f.write('"%s"' %(options['Title'])  +  tab + str(stacks.get_numElements()) + tab + str(len(self)) + le)
        
            #f.write(ts[9])  # comment: PlotType...
            f.write("PlotType (0-5); Plot Depths: Xmin, Xmax(Ang.) [=0 0 for Viewing Full Target]"+le)
            f.write(tab + str(options['PlotType']) + tab + str(options['PlotExtents'][0]) + tab + str(options['PlotExtents'][1]) + tab + le)
        
            ## List of Elements:
            #f.write(ts[11])  # comment: Target Elements...
            f.write("Target Elements:    Z   Mass(amu)"+le)
            for ii,e in enumerate(stacks.elements):
                f.write("Atom %i = %s =      %f  %f" %(stacks.elidnum[ii], e.name, e.elnum, e.mass)    + le)
        
            ## List of Target layers:
            #f.write(ts[13])  # comment: Layer...
            f.write("Layer   Layer Name /               Width Density    ----Elements go here----"+le)
            #f.write(ts[14])  # comment: Numb.  ...
            f.write("Numb.   Description                (Ang) (g/cm3)    Stoich  Stoich  Stoich  Stoich  Stoich"+le)
            # 1      "GaSb"           4830  6.294      .5      .5       0       0       0
            prevelt = 0 # position of previous element
            for i,l in enumerate(stacks.stack):
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
            
                for ii in range(stacks.get_numElements() - numelts ):
                    # fill rest of elements with 0 mole frac.
                    fstr += '0' + tab
                
                fstr += le    
                f.write(fstr)
            #end for(Layers)
        
            #f.write(ts[16])  # comment: Target layer phases...
            f.write("0  Target layer phases (0=Solid, 1=Gas)"+le)
            f.write(' ')
            for l in stacks.stack:
                B = 1   if  l.isGas else  0
                f.write( str(B) + ' ' )
            f.write( le )
        
            #f.write(ts[18])  # comment: Target compound corrections...
            f.write("Target Compound Corrections (Bragg)"+le)
            f.write(' ')
            for l in stacks.stack:
                f.write( str(l.compoundCorrection) + '   ' )
            f.write( le )
        
            #f.write(ts[20])  # comment: atom displacements...
            f.write("Individual target atom displacement energies (eV)"+le)
            for e in stacks.elements:
                f.write( tab + str(e.displacement)  )
            f.write( le )
        
            #f.write(ts[22])  # comment: atom lattice binding...
            f.write("Individual target atom lattice binding energies (eV)"+le)
            for e in stacks.elements:
                f.write( tab + str(e.binding)  )
            f.write( le )
        
            #f.write(ts[24])  # comment: atom surface binding...
            f.write("Individual target atom surface binding energies (eV)"+le)
            for e in stacks.elements:
                f.write( tab + str(e.surfbinding)  )
            f.write( le )
        
            #f.write(ts[26])  # comment: Stopping power verions
            #f.write(ts[27])  # Stopping power verions
            f.write("Stopping Power Version (1=2011, 0=2011)"+le)
            f.write(" 0 " +le)
        
            f.close()
        #end for(allstacks)
    #end output()
    
    
    def splitstack(self, numlayers=30):
        '''Split up this stack into multiple smaller stacks, according to the maximum number of layers specified by _SRIMMaxLayers.
        
        Parameters
        ----------
        numlayers : int, optional
            How many layers should the resulting lists contain?  Defaults to 30.
        
        
        Returns
        -------
        
        ListOfStacks : list
            Returns a List conatining the split-up sub-stacks, where each sub-stack has no more layers than specified by `numlayers`.  The sub-stacks can be output as normal.
            This function is used by Stack.output to automatically split-up a Stack with more layers than SRIM can handle, and create multiple output files accordingly.
            
        Examples
        --------
        
        For example, if the original Stack was created via:
        
            >>> Orig = Stack( 18*( Mat1(thickness) ) )  # 18 layers
            
        And it was split with 
        
            >>> NewStacks = Orig.splitstack(numlayers=5)    # want 5 layers max
            
        The resulting list would look like so:
        
            >>> NewStacks
            :   [StackWith5Layers, StackWith5Layers, StackWith5Layers, StackWith3Layers]
        
        '''
        if DEBUG(): print "numlayers=", numlayers
        n = max(1, int(numlayers))   # sanitize numlayers, to minimum of 1
        newstacks =  [self.stack[i:i + n] for i in range(0, len(self.stack), n)]
        # newstacks is list of lists of Layers:
        #   [[<pyTRIMSetup.Layer object>, <pyTRIMSetup.Layer object>], [<pyTRIMSetup.Layer object>]....]
        return  [Stack(x)   for x in newstacks]  # converts each list element into a Stack object
        
    #end splitstack()
    
#end class(Stack)
    
        


