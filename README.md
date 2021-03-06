
pyTRIMSetup
===========

A Python module for generating a TRIM.IN input file for running TRIM.exe by James Zeigler (see [srim.org](http://srim.org)).
SRIM/TRIM is an ion-implantation monte-carlo simulator.

After generating an output file with your python script, save the file into the SRIM/TRIM directory, and rename it to `trim.in`.

You can run your simulation in two ways.  The preferred method is to Run `SRIM.exe`, and click `TRIM Calculation`/`Restore last TRIM data`, which loads the `TRIM.IN` file into the SRIM GUI.  You can then click `Save Input and Run TRIM`, which will subject your input to the input checks, such as "merging identical elements" (useful because TRIM can fail with a run-time error if too many elements are present, and merging elements dramatically speeds up the simulation).

Alternatively, you can run `TRIM.exe` (skipping the setup GUI), which will take all settings directly from your `TRIM.IN` file and immediately start the simulation.  TRIM will load the settings from your generated file and start the simulation.
This bypasses the SRIM.exe graphical interface entirely.

Plan to develop the automatic launch of TRIM.exe with the generated TRIM.IN file & stitching together many-layered targets via the TRANSMIT.txt output option.

####Contact

Feel free to add issues/feature requests, or even better, Fork the `git` repository and submit your own updates (see this [how-to](http://kbroman.org/github_tutorial/pages/fork.html))!

Nov. 2015, [Demis D. John](mailto:demis.john@gmail.com)

---------------------------------------------------------------
## Example:

###Import the pyTRIMSetup module

    from pyTRIMSetup import *   



###Setup target materials

A Material object is created as so:

  *NewMat = Material(  [List_Of_Element_Abbreviations],  [List_Of_Mole_Ratios],  Density)*

    GaSb = Material(['Ga','Sb'], [0.5,0.5],   3.657)
    AlAs = Material(['Al','Ga','As'], [0.3,0.3,0.4],   2.758)
    AlSb = Material(['Al','Sb'], [0.5,0.5],   1.97)



###Setup Target layer Stack

  The Stack's Layers are created from top-to-bottom, by adding materials together, while providing a thickness (*in Angstroms*) for each layer.

    target = Stack(  GaSb(1500) + AlAs(750) + GaSb(2000) + AlSb(2500) )     # top to bottom

Since the underlying structure of these commands uses python Lists, we can do some interesting operations with them, like so:

    # One period of AlAs & GaAs is multipled to repeat it:
    repeatingpart = 3 * (  AlGaAs(150) + GaAs(110)  )   
    
    # Insert this into the main Stack:  
    target = Stack(  GaAs(1500) + repeatingpart + AlAs(2500)  )     # top to bottom



###Setup ion to implant

The Ion to implant is defined as

  *Ion(  Element_Abbrv, Energy_keV, Angle_degrees )*
  
    target.implant(    Ion('H', 10, 7)    )     # Ion(ElementAbbrv, Energy_keV, Angle_degrees)



###Setup a dictionary of simulator options:

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



###Generate the output TRIM.IN file

Output the file at specified path/filename, passing the above options dictionary.

    target.output('TestOutput.in', options=options, overwrite=True)


The file 'TestOutput.in' can now be copied to TRIM.IN in the SRIM/TRIM directory.
If you run TRIM.exe, it will pick up all the simulation settings from this file.

Alternatively, run SRIM.exe, and choose "TRIM Calculation" & "Restore Last TRIM Data", which will load your simulation into the GUI.  This is useful because then, when you choose "Run TRIM", it will subject your layers to the input checks, the most useful of which is "Merge Identical Elements".  It is recommended you do this, as otherwise TRIM can crash with a run-time error when too many elements are present (it also runs horribly slowly).

A planned future feature is to automatically copy the output file to TRIM.IN & run TRIM.exe for you - optionally stitching together long simulations by running multiple TRIM's.

