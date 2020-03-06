========
CONTENTS
========

	1 - RUNNING THE PROGRAM
	2 - INPUTS
		poskms.txt
		numstepsFAST.txt & numstepsSLOW.txt
		vmd.xyz
	3 - OUTPUTS
		Graphs
			Solar System simulation: total energy vs time
			Moon-Earth Separation: separation vs time
		Files
			timeperiodfile.txt
			apoperifile.txt

===================
RUNNING THE PROGRAM
===================

Instructions on how to run the program:
In bash the files needed to be loaded so that the program runs are, in that order:
		VelocityVerlet3D poskms.txt numstepsFAST.txt vmd.xyz 

The units used in the simulation are km, kg and days

==================
      INPUTS
==================

poskms.txt
================
Holds initial positions and velocities of each paticle, with a unique label.

Format
------
label x y z vx vy vz
Type: label - int || x,y,z,vx,vy,vz - float
Units: km & km/day


numstepsFAST.txt & numstepSLOW.txt
================
Holds simulation parameters, including G, time step, and simulation length in steps.

Format
------
G dt time_end
Type: G - float || dt, time_end - int (*)
Units: G - [ km^3/(kg*days^2)] || dt, time_end - days 

*Can also be inputted as other types, but will be converted to ints within the code

NOTE: there is an option to run numstepsSLOW.txt or numstepsFAST.txt. 
FAST will use dt - 4, time_end - 200000 (approx 2 pluto orbits)
SLOW will use dt - 1, time_end - 100000 (approx 1 pluto orbit)

SLOW will make the simulation take longer to finish but will produce more accurate results, specifically when it comes to the revolution of the moon around the Earth. CAUTION: It will also make the VMD file significantly larger (20MB -> 41MB).

vmd.xyz
================
Name of the file in which your VMD data will be output to.
Can be renamed to whatever suits the user. CAUTION: ensure it ends in .xyz or else it will break!

Format
------
As this is an empty/non-existent file on input, so this has no format.
After the VMD data is written to this file, it uses the format required by VMD
described in the CompMod Planetary Simulation instruction manual.

NOTE: This file size can vary largely based on your numstepsFAST/SLOW.txt file, so be careful when fiddling with the simulation parameters.

==================
     OUTPUTS
==================

Graphs
===============
a) Solar System simulation: total energy vs time

Plot of the energy of the system, against the time of the simulation. 
Units are in 1.3996E-4 Joules and days.

b) Moon-Earth Separation: separation vs time

Plot of the separation of the Moon & Earth, against the time of the simulation. 
Units are in km and days.

Files
===============
a) time_periodfile.txt 

Contains the periods of all the planets and Halley's comet around the Sun.
Units used: Earth years; This is to make comparison to literature data easier but also it is a more natural way of portraying time scale instead of Earth days.

b)apoperifile.txt 

Contains the peri and apoapsis of all planets and Halley's comet around the Sun.
Units used: Km

c) energy.txt

Contains the total energy of the system at each timestep in a list
Units used: 1.3996E-4 Joules and days.

***Note: Both apoperifile and time_periodfile have additional accounts to the time period and apo- and periapsis of the moon around EARTH. 
