************************
# tracer2D Manual
************************

The purpose of this application is to calculate the paths of mass-less tracer particles started from positions defined in a text file. The input flow field should be periodic, since when the tracing time surpasses the time length of the input flow field it is reused periodically from the beginning.


## How to set it up?

tracer2D uses the coordinate system of the geometry image, thus, every starting position or outlet definition should use pixel coordinates ([0,0] is in the upper left corner).
The main setup is contained in an .ini file, as ususal.
Additional necessary files are particle starting points (containing two columns of pixel coordinates: x y),
and outlet definitions again in a text file containing four values per line, easch corresponding to an outlet where particles can leave.
Note: outlets are defined as a rectangular regions with two corner points as: x1 y1 x2 y2 where x2 > x1 and y2 > y1. 

The sample tracing data is set up to work with the results of the 'art_aneur' simulation.


## How to run it?

### Windows

1.) Check medFlow2D/data_tracer/trace_setup.ini for a sample setup.
2.) Currently using the command-line, from medflow2D/bin:
	tracer2D.exe ../data_tracer/trace_setup.ini
