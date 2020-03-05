Instructions on how to run the program:
In bash the files needed to be loaded so that the program runs are, in that order:
		VelocityVerlet3D poskms.txt numstepsFAST.txt vmd.xyz 

Note: there is an option to run numstepsSLOW.txt instead. This file will make the simulation take longer to finish but will produce more accurate results, specifically when it comes to the revolution of the moon around the Earth.

The units used in the simulation are km, kg and days
For G =  [ km^3/(kg*days^2)]

The output of the simulation is:
a) Plot of the energy of the system, against the time of the simulation. 
b) time_periodfile.txt --> Contains the periods of all the planets and Halley's comet around the Sun.
Units used: Earth years; This is to make comparison to literature data easier but also it is a more natural way of portraying time scale instead of Earth days.
c)apoperifile.txt --> Contains the peri and apoapsis of all planets and Halley's comet around the Sun
Units used: Km

***Note: Both apoperifile and time_periodfile have additional acoounts to the time period and apo and periapsis of the moon around EARTH. 