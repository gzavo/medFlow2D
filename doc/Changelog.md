## 2018.03.03
- updated most python scripts to python 3 compatibility
- added small karman vortex street example

## 2015.08.31
- tweaks to fractaldimGUI.py

- - - - - - - - - - - Version 1.2 - - - - - - - - - - -

## 2015.08.30
- tweaks to fractaldimGUI.py
- tweaks to positionPlot scripts
- win64 / intel maci64 / linux64 binaries for both medFlow2D and tracer2D

## 2015.08.29
- fixed linux compatibility

## 2015.08.28
- fractaldimGUI.py tweaks
- tracer position plot tweaks

## 2015.08.23
- Fixes in coloring mode
- Fixes in water ageing dynamics (wrong PSRhoMean)
- Small fixes in Monitor.py for accurate plotting of water ageing
- Improvements in the theory guide (e.g., water aging...)
- Minor directory restructuring
- Added plot tool for particle position visualization

## 2015.08.01
- Small changes in the _Manual_.
- GUI python code for fractal dimension analysis.

- - - - - - - - - - - Version 1.1 - - - - - - - - - - -

## 2015.07.30
- Implemented preliminary fluid ageing algorithm (uses the coupled field)
- BREAKING CHANGE: Changed the color scaling routine to match ageing data too: imgRhoScale and imgPSScale now corresponds to the tota width of data scale. E.g., for density is typically in [0.95; 1.05], imgRhoScale = 0.1 will give a good scale.

## 2015.07.28
- Small fixes to the chaotic analysis package (only the plotting part was modified)
- Made tracer2D OpenMP parallel (thread parallelism implemented too)

## 2015.07.24
- Small fixes in the documentation and the validations

## 2015.07.23
- Some file restructure in tools
- Added genGeom_straight.py, a straight geometry generator
- Updated monitor.py (log scale plotting)
- Fixed several bugs in monitor.py
- Added new status at contributors: "developer" (see: Credits.md)

## 2015.07.17
- Streamlining the documentation a little
- Updated _Theory Guide_
- Fixed ETA computation to remain stable
- New simulation template -> artAneur_curved (generated with genGeom.py)

## 2015.07.16
- Changed main icon file
- added genGeom.py to tools for generating parametric aneurysmal geometry

- - - - - - - - - - - Version 1.0 - - - - - - - - - - -

## 2015.07.13
- Fixed tracer2D thread race condition bug
- Small updates to the GUI
- Since all major parts are functional: bump version to 1.0

## 2015.07.08
- Small fixes to status messages
- Fixed setting for Smagorinsky constant
- Implemented forced MRT dynamics (e.g., MRT for porous material)
- Extended the TheoryGuide

## 2015.07.07
- Changed MRT relaxation parameters for higher stability
- Small tweaks to GUI behaviour
- Small tweaks to command file handling (e.g., nicer cleanup on user termination)

## 2015.07.06
- Upgraded to build with GCC 5.1
- Reimplemented MRT method in a parallel manner (it can now work with multiple workers)
- plotData.py extended with figure title and captions
- Changed from OpenMP to thread parallelization to improve compatibility across devices
- Several tweaks for linux compatibility

## 2015.07.03
- Added possibility to set dx resolution of the input files (du/dy is dependent on that information to give SI results).
- Small tweaks to the documentation.

## 2015.07.02
- Fixed BUG at saving data files but not the images

## 2015.07.01
- Changes y-axis direction of contour plots
- added tracer2D program (to be integrated to the GUI later?)

## 2015.06.25
- Added validation setups and sample results to 'data/validations'

## 2015.06.22
- plotData.py is merged with the main control GUI (e.g., the 'Post' tab of the interface)

## 2015.06.21
- Implemented data compression routines (using miniz library)
- Implemented PNG image compression routines (using miniz library)
- All image outputs are now in png format (smaller)
- Data output is now zipped by default
- monitor.py and plotData.py have been improved to support zip compression

## 2015.06.20
- Implemented logging system
- The simulation saves the setup.ini to the output folder at the beginning of the simulation

## 2015.06.19
- Restructured output file system
- BREAKING CHANGE: setup.ini tags corresponding to file output has been simplified
- Updated all the setup.ini files

## 2015.06.16
- Implemented arbitrary inlet velocity directions (needed for lid driven cavity)
- Implemented a new switch to choose inlet profile
- Moved monitor.py to gui directory
- Small fixes to monitor.py
- Restructured setup.ini (karman example only so far)
- version bumped to 0.8

## 2015.06.15
- Implemented Smagorinsky turbulence model for the BGK dynamics
- Added setup.ini keywords 'useSmagorinsky' and 'smagorinskyC' (see karman sample)
- Small fixes and improvements in the editor module.
- Improvements in the 'TheoryGuide'

## 2015.06.13
- Reimplemented the regularized boundary conditions (they are less accurate, but stable up to 10x higher Re)
- Updated the TheoryGuide heavily
- New settings under [Parameters]: useRegularized to enable the new boundary condition
- small fix to monitor.py to make sure the input text file is closed even on error
- Improved the setup editor a little bit.

## 2015.06.11
- Reimplemented the BGK dynamics, the Zou-He boundary conditions, the propagation method
- version bumped to 0.7

## 2015.06.10
- Started new documentation method
- Added a theory guide (a markdown editor is advised. e.g., haroopad)
- Started reworking the coordinate systems and the boundary conditions

## 2015.05.22
- New tool to visualise maximal Pcoag values during the simulated time
- Small fixes in the other python tools

## 2015.04.29
- Changed the naming convention in the coag_*.txt outputs.
- Fixed a few glitches in coagPlot.py.
- new simulation preset: data/pipe2

## 2015.04.28
- BREAKING CHANGE: Reynolds number calculation changed...again...: v_avg = 3/2 * v_max. Re=v_avg*L/nu.
- Added plotCoag.py to tools. (To plot coagulation probabilities.)
- Added plotVector.py to tools. (To plot velocity and margination force.)

## 2015.04.13
- monitor.py is extended with autoscrollbars
- new simulation preset: data/aneur_coag

## 2015.04.09
- The program no longer saves the density image (it is useless anyways)
- Added more (error) messages around the ini file loading
- BREAKNG CHANGE: renamed imgSaveInterval to saveInterval as it controls all saving frequency not just the image saving

## 2015.03.26
- Detecting OpenMP support for the used compiler

## 2015.03.23
- Added glycocalyxWidth to ini file setup to account for the possible lack of near wall margination forces (RBCs typically wont go that close to the walls)
- Major code refactoring to prepare for double buffered computation and fused collide and stream dynamics

## 2015.03.17
- Refactored coagulation code to use proper names
- BREAKING CHANGES: input ini file renamed parameters under [thrombus]:
	coagProbDecr -> rhoADPDecr
	defCoagProb -> defRhoADP
- new parameter in the ini files under [thrombus]: coagCoeff

## 2015.03.12
- Fixed BUG that caused program crash when changing from warmup to simulation
- Protected fluid sites (cyanide blue) now don't experience margination forces (thus reducing the inlet/outlet error effects)
- Made simulations more efficient by eliminating unused wall cells
- Corrected image output colours
- Updated CMAKE file for new cmake versions

## 2015.03.10
- Refined inner calculation of dimensionless reference system
- Corrected pressure output values in [Pa]
- Corrected porous constant calculations

## 2015.03.06
- BREAKING CHANGE: Changed Reynolds number to be calculated with the average flow velocity across the inlet instead of the maximal. (With parabolic profile U_avg = U_max/2.)

## 2015.03.05
- Added error handling for the transient boundary data file
- Corrected small bugs around the .command file
- Added basic GUI to control the simulation
- Added primitive text editor to GUI
- Added timer to the GUI to check if subprocess is alive

## 2015.03.04
- Added possibility to create protected fluid cells (new colour)
- Added save option for coagulation probabilities
- Added info to 'tools' on how to plot the coagulation probabilities
- Added new example (artaneur_coag)
- The manual is extended with additional info

## 2015.03.03
- Added data saving capability to plotData.py

## 2015.03.02
- Fixed plot aspect ratio in plotData.py
- Fixed typo in plotData.py

## 2015.03.01
- Bumped version to 0.6
- Several small bugfixes in error handling and in messages
- Switched to parallel execution (OpenMP)
- Print performance in MLUPS (Million Lattice Updates Per Second)
- Estimating remaining time in the given simulation phase
- Added post processing script templates to 'tools'
- Updated all ini files to be correct
- Added new simulation (karman)
- Fixed aspect ratio in the plots of monitor.py

## 2015.02.24
- Fixed the BUG responsible for the huge memory consumption
- Changed monitor.py to tabbed view
