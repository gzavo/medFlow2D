# Manual for medFlow2D

-----------

## Welcome

### About
medFlow2D is a CFD software capable of solving various two-dimensional flows. It was primarily designed for educational and scientific use, thus easy setup of the flow problems is an expected behaviour.
The solver is built on the lattice Boltzmann method (see _Theory Guide_ for deatils).

### Features
- Dynamics: BGK, MRT
- Smagorinsky turbulence model
- Boundaries: Zou-He, Regularized
- Macroscopic porous material modeling
- Non-Newtonian blood model
- Coupled passive scalar field
- Macroscopic margination modeling
- Thrombus formation model
- Platelet margination model
- Fluid age calculation
- Particle path tracing
- Fractal analysis
- Parallel execution
- Complete toolsuite (Pre/Solver/Monitor/Post)
- Cross-platform


-----------

# Quickstart for Windows

Tested under Windows 7/8.1.

## Installation

Requirements, solver: the solver has no external dependency.
Requirements, GUI: needs a working python installation with numpy, matplotlib, ImageTk (e.g., Pillow)

Recommended installation:

- Install Anaconda Python distribution ( https://store.continuum.io/cshop/anaconda/ )
- Install Pillow:
    > conda install pillow

- Start solver control with 'startGUI.bat'
- Check the bundled cases in the _data_ folder

-----------
# Quickstart for Linux

Tested on Ubuntu 14.10.

## Build under Linux

Extract medFlow2D sources, then:

> cd medFlow2D
> mkdir build
> cd build
> cmake ../ -DCMAKE_BUILD_TYPE=Release
> make

## Start solver control

> sh startGUI.sh

- Check the bundled cases in the _data_ folder

-----------


## Set up a new simulation

Use ini definition files to set up the simulation (see the bundled samples).
The geometry is described by a ppm "paintable" bitmap file (for the colour coding see sections below).

### Color coding for the geometry

R=0   G=0   B=0   -> Fluid
R=0   G=255 B=255 -> Protected fluid (no coagulation or other special process allowed here)
R=255 G=255 B=255 -> Wall
R=255 G=0   B=0   -> Inlet (velocity boundary, one one of these are allowed)
R=0   G=255 B=0   -> Outlet (pressure boundary)
R=0   G=0   B=255 -> Porous material

### Tips for editing the geometry

- The simulation accepts PPM (Portable PixelMap) format. Several free image viewer can convert from and to this format, e.g., http://www.irfanview.com/ .
- Geometry can be edited in any image editor, but I would recommend Paint.Net ( http://www.getpaint.net/index.html ) under windows and GIMP ( http://www.gimp.org/downloads/ ) on any other operating system.
- Be sure to turn off anti-aliasing in the chosen editor (otherwise unwanted shaded colors will appear in the image)!

### Available boundary conditions

- The simulation can have one inlet opening, and arbitrary number of outlet openings.
- The value of the inlet scale-function should scale between [0; 1]
- If the scale-function is shorter in time than the total simulation time, the program will start the scale function from the beginning (periodically).
- The values between two scale-function value are interpolated linearly
- Inlet nodes can only be located on a single side, without encorporating corner nodes!

- only one inlet is allowed (!!!)
    - the inlet is defined as velocity inlet with an automatically generated parabolic profile
- arbitrary number of outlets are allowed [as constant pressure BCs.]

------------

## Install required softwares

The solver itself has no external dependency.
Since it is a command line interface application, you can use it through the command line.

The software package provides further basic tools for visualization, simulation management and monitoring.
However, these tools require the presence of a python interpreter, and a few python packages (tk, pillow).
The recommended python interpreter can be downloaded from:
Python - https://store.continuum.io/cshop/anaconda/
During the install select "Add to system PATH".
You can check if the setup has been completed correctly by asking for the version of the compiler by executing the following in a command line:

python --version

After Anaconda setup is completed, open a new terminal to install the only missing package, Pillow.
To do this, enter the following command to the command line:

pip install Pillow
(You might need a new version of Pillow -> pip install Pillow --upgrade)


## Run and control a simulation using the Solver manager

If all required softwares are present on the system you can simply execute 'startGUI.bat' (Windows) or 'startGUI.sh' (Linux/Mac) to start the solver manager application.



## Alternatively you can run and control the simulation from command line

This can be useful if you execute the simulation and the monitoring tools on different machines.

1. Extract the program somewhere.
2. Pick a simulation in folder 'data'.
3. Edit the .ini file for that simulation (e.g., 'setup_stent.ini').
4. Make sure the two output paths are correct and exist in your system:
'base_name' under [output] will be the directory for the simulation output,
and 'currentInfo' under [monitoring] will hold the monitor files.
3. In a command line navigate to the 'bin' directory.
5. Run the simulation by executing:
    a. Windows: medFlow2D.exe <setup.ini>
    b. Linux/Mac: ./medFlow2D <setup.ini>
    where <setup.ini> contains the description for the simulation.
6. Start the monitoring software:
    a. Open a new command line.
    b. Navigate to the 'tools' directory.
    c. Execute: python monitor.py <currentInfo path you jast set>.
7. You can stop the simulation by pressing CTRL+c in the command line where the simulation runs.

## Compilation on Windows

The solver is written in ANSI C. It can be compiled with any standard compatible compiler. The recommended compiler can be downloaded from:
GCC - http://tdm-gcc.tdragon.net/
The 64bit version is most likely to be the one you want to download.
During the install select "Add to system PATH".
You can check if the setup has been completed correctly by asking for the version of the compiler by executing the following in a command line:
gcc --version
The compilation directives are stored in a CMake file. Thus, you can compile the source using CMake. ( CMake - http://www.cmake.org/ ). Or you can use an IDE that has CMake embedded in. I recommend CLion:
C++ IDE - https://confluence.jetbrains.com/display/CLION/Early+Access+Program

### Windows FAQ

Q: How to access the command line in windows?
A1: Start menu -> Run (or Win+r keys) then execute 'cmd'
A2: If you use Total Commander navigate to the desired directory, and execute 'cmd' there.



## Compilation on Linux

Compilation requires a recent version of the GNU compiler suite.

In the root directory of medFlow2D:

> cd bin
> cmake ../
> make
> ./medFlow2D <some_ini_file>



## Additional general info
- There can be only one inlet and arbitrary number of outlets
- Re is calulated with the average of the inlet velocty and the inlet width



## Advanced informations

For deatils consult the _TheoryGuide_.

### Coagulation

- Fluid cells next to a wall are defined with bgkPS instead of bgkForcedPS to increase warmup stability.
- effect of glycocalyx can also be modelled using layers of bgkPS instead of bgkForcedPS


------------

Known issues
============

- [FIXED] Porous material only works for BGK mode (develop alternative for MRT)
- [FIXED] Instable platelet density
- [FIXED] Check upper pressure boundary condition (seems incorrect)
- [FIXED] Monstrous memory consumption
- [FIXED] Scale ADP [kg/m^3] with lattice size
- [FIXED] Something is wrong with the Zou-He boundaries, check the implementation (velocity drop at inlet, probably due to pressure increase????)

Development TODO
================

- Wall normal and WSS calculation
- Test double buffered dynamics for performance
- Allow for a constant translation of platelet concentration
- Implement factorized central moments
- Compute proper inlet flow based on inverse Womersley
- optimize outlet BC (zero gradient output/windkassel model approx.?)
- Driest damping for Smagorinsky turbulence model
- Improved error handling (for file operations, wrong parameters, etc...)
