
== pyTRIMSetup ==

A Python module for generating a TRIM.IN file for running TRIM.exe by James Zeigler (see srim.org).
SRIM/TRIM is an ion-implantation monte-carlo simulator.

After generating an output file with this script, save the file into the SRIM/TRIM directory, rename it to TRIM.IN and run TRIM.EXE.  It will load the settings from your generated file and start the simulation.
This bypasses the SRIM.exe graphical interface, and allows for repeating loops of target layers etc.

Nov. 2015, Demis D. John

---------------------------------------------------------------
Example:

*import the pyTRIMSetup module*

    from pyTRIMSetup import *   


*Setup target materials*
A Material object is created as so:
  NewMat = Material(  [ListOfElements],  [ListOfMoleRatios],  Density)

    GaSb = Material(['Ga','Sb'], [0.5,0.5],   3.657)
    AlAs = Material(['Al','Ga','As'], [0.3,0.3,0.4],   2.758)
    AlSb = Material(['Al','Sb'], [0.5,0.5],   1.97)


*Setup Target layer Stack*
  The Stack's Layers are created from top-to-bottom, by adding materials together, 
  while providing a thickness (Angstroms) for each layer.

    target = Stack(  GaSb(1500) + AlAs(750) + GaSb(2000) + AlSb(2500) )     # top to bottom


*Setup ion to implant*
The Ion is defined as
  Ion(  ElementAbbrv, Energy_keV, Angle_degrees )
  
    target.implant(    Ion('H', 10, 7)    )     # Ion(ElementAbbrv, Energy_keV, Angle_degrees)


*Setup a dictionary of simulator options:*
  These will eventually be made as default options, but right now you need to specify each one.

    options = {}
    options['Title'] = 'Testing the script' 
    options['NumIons'] = 5000       # Max number of ions to simulate - simulator stops at this many ions
    options['AutoSaveNum'] = options['NumIons']     # Save the simulation at this interval
    options['SimType'] = 2      # 1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons.  See SRIM Menu for more info.
    options['RandomSeed'] = 0
    options['Reminders'] = 0
    options['DiskFiles'] = [0,0,0,0,0,0]    # booleans for: Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file
    options['PlotType'] = 5     # 5: no plot.  See the SRIM Plots menu.
    options['PlotExtents'] = (0,0)      # Xmin, Xmax = 0,0 for automatic



*Generate the output TRIM.IN file* with specified path/filename, passing the above options dictionary.

    target.output('TestOutput.in', options=options, overwrite=True)

    '''
    The generated file can now be copied to TRIM.IN in the SRIM/TRIM directory.
    If you just run TRIM.exe, it will pick up all the simulation settings from this file.
    '''
