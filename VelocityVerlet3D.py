"""
Author: Panayiotis Panayiotou, John Waiton
StudentID:s1715899, s1739002
Project: N-body Simulation of the Solar System
"""
import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D

# add the Gravitational constant
# Units used: km, kg and days
# For G, in units km^3/(kg*days^2)

## files neeeding loaded in, particle positions, velocities & constants (numstep, dt, G)

## Particle positions taken from 28th February 2020 from NASA's site via CompMod instructions

def fileread(arg):
    """
    Method to collect values from input files provided
    and register correct files to read/write from/to.

    :param arg: list of all input and output files;
                initial positions & velocities of particles, system constants
                (dt, length of simulation, gravitational constant)
    """
    if len(arg) != 4:
        print("Wrong number of arguments.")
        print("Usage: " + arg[0] + " " + arg[1] + " " + arg[2] + arg[3] + " " + " Another file")
        return
    else:
        # Collect input filenames
        inputfile_name1 = arg[1]
        inputfile_name2 = arg[2]
        # register vmd file name
        vmdfile = arg[3]

    # Open input files
    positions = open(inputfile_name1, "r")
    int_param = open(inputfile_name2, "r")

    # Collect simulation parameters from 2nd input file
    file_handle2 = int_param
    G,dt,time_end = Particle3D.from_file2(file_handle2)
    # Convert from string to float and int to ensure run time is divisible by dt
    G = float(G)
    dt = int(dt)
    time_end = int(time_end)
    return G,dt,time_end,vmdfile,positions


def newton(p1, p2, m1, m2):
    """
    Method to return potential energy and force
    of particle using Newtonian mechanics

    :param p1 and p2: Particle3D instances
    :param m1 and m2: masses of the objects under gravitational potential
    :return: potential energy of particle as float & force acting on particle as Numpy array
    """
    sep_scalar_vec = Particle3D.seperation(p1,p2)
    sep_scalar= np.linalg.norm(sep_scalar_vec)
    potential = -(G*m1*m2)/sep_scalar
    force = -(G*m1*m2*sep_scalar_vec)/(sep_scalar**3)
    return force, potential


# A function is created to condense the code as this block of code is repeated throughout.
def baseloop(particle):
    """
    Method to process any function across
    all planets and apply results to matrix

    :param particle: Particle list
    :return: Matrix of applied force/energy in suitably sized matrix
    """
    # example of use: totalF = baseloop(particle,force_Newton)
    # Initial conditions

    # Creates w lists, each with h elements. Effectively,creating a list of lists
    # For our use, each list is a particle, each element is a force
    # to initialise our matrices we use list comprehension.
    w, h = (len(particle)),(len(particle)-1);
    # fmatrix will hold the forces, pmatrix holds the potentials.
    fmatrix=[[0 for x in range(w)] for y in range(h)]
    pmatrix=[[0 for x in range(w)] for y in range(h)]
    # begin basic looping system, allows for looping across all planets
    for f in range(0,len(particle)-1):
        # No particle, break loop. Failsafe
        if particle[f] == None:
            break
        else:
             for g in range(0,len(particle)-1):
                # If we are not testing against ourselves, apply function
                 if f != g:
                     if particle[g] == None:
                         break
                # To stop us calculating already existent results
                     if f > g:
                         fmatrix[f][g], pmatrix[f][g] = -fmatrix[g][f], -pmatrix[g][f]
                     else:
                         fmatrix[f][g], pmatrix[f][g] = newton(particle[f].position,particle[g].position, particle[f].mass,particle[g].mass)

    # Apply these to find the total forces and potentials of each list (planet)
    ftotal = [sum(x) for x in fmatrix]
    ptotal = [sum(x) for x in pmatrix]
    return ftotal, ptotal


def apoperi(particle, apo, peri):
    """
    Method to update the apo and periapsis of each particle
    with respect to the sun for all particles except the moon
    which is with respect to Earth

    :param particle: Particle list
    :param apo: Existing apoapses list
    :param peri: Existing periapses list
    :return: Lists of all particles' updated apo- and periapses
    """
    # Loop for all planets except sun
    for f in range(1,len(particle)-1):
        # All separations excluding the moon. Special consideration for the moon.
        if f != 4:
            # Calculate separation
            sep = np.linalg.norm(Particle3D.seperation(particle[0].position,particle[f].position))
            # Determine if new apoapsis
            if apo[f] == None:
                apo[f] = sep
            elif (sep > apo[f]):
                apo[f] = sep

            # Determine if new periapsis
            # If no periapsis exists
            if peri[f] == None:
                peri[f] = sep
            elif (sep < peri[f]):
                peri[f] = sep
        # For the moon, relative to the Earth
        else:
            # Calculate separation, assuming Earth is particle 3
            sep = np.linalg.norm(Particle3D.seperation(particle[3].position,particle[f].position))
            # Determine if new apoapsis
            if apo[f] == None:
                apo[f] = sep
            elif (sep > apo[f]):
                apo[f] = sep

            # Determine if new periapsis
            if peri[f] == None:
                peri[f] = sep
            elif (sep < peri[f]):
                peri[f] = sep
    return apo, peri

def dot_prod(a,b):
    """
    Method to calculate the dot product of two vectors
    and normalise said result.

    :param a: Vector 1 in our dot product
    :param b: Vector 2 in our dot product
    """
    # Calculate dot product and normalise
    dot_prod = np.dot(a,b)
    value = dot_prod/(np.linalg.norm(a)*np.linalg.norm(b))
    return value


########################## BEGIN MAIN CODE ########################################
def main():
    # Files are read into the system before main(); as can be seen at the end of this code.
    # Two inputfiles were chosen, one that has information about the particles
    # while the other has information about the interaction and the amount of iterations carried.
    # One output file is chosen by the user, will write the vmd to said file (ensure it ends in .xyz!)
    # code written; python3 Velocityverlet.py poskms.txt numsteps.txt vmd.xyz

    ###################### PARTICLE CREATION #####################

    # Initialise particle list
    particle= []

    # Extract information from file about the positions of the particles and their velocity, create Particle3D instances stored in list
    file_handle1 = positions
    i=1

    particle.append(Particle3D.from_file1(file_handle1))
    # create loop to allow formation of each particle.
    while True:
        particle.append(Particle3D.from_file1(file_handle1))
        # when particle is empty (list reading finished) end while loop
        if particle[i] == None:
            break
        i+=1
    # show particle number and completion
    print("Particle Processing Complete, particles processed: " + str(len(particle)-1))

    ##################### VARIABLE INITIALISATION ###########################

    ## This section doesn't include all variable initialisation, but ones that are used in time integration extensively
    # Initialise time.
    time = 0.0

    # Initial force & energy
    # Use the numpy module to create a matrix that will store the values of the force
    totalforce=[]
    totalforcenew = []
    totalenergy = []
    # Initialise data lists for plotting later
    # Lists that will contain the time evolution of the position of both particles and their seperation as well as their energy will be taken
    # Since the position of particles is a numpy array, we need to find their norm so we can plot it
    time_list = [time]
    # Used for plotting the moon-earth separation
    moonsep_list = []
    # Labelling each particle, separate from particle.label due to discrepancy in value (0->1)
    plabel = []                                 #Particle No.           #Read in from position file as
    plabel.append("Sol")                        #particle 0             1
    plabel.append("Mercury")                    #particle 1             2
    plabel.append("Venus")                      #particle 2             3
    plabel.append("Earth")                      #particle 3             4
    plabel.append("Moon")                       #particle 4             5
    plabel.append("Mars")                       #particle 5             6
    plabel.append("Jupiter")                    #particle 6             7
    plabel.append("Saturn")                     #particle 7             8
    plabel.append("Uranus")                     #particle 8             9
    plabel.append("Neptune")                    #particle 9             10
    plabel.append("Pluto")                      #particle 10            11
    plabel.append("Halley's Comet")             #particle 11            12


    # Calculate initial total force & total potential energy
    totalforce, potenergy = baseloop(particle)

    ##################### WRITING TO FILES #####################
    # Open Files
    efile = open("energy.txt","w")
    partfile = open(vmdfile,"w")

    # Write total energy to Files
    # Resetting energy value
    energy = []
    # Write break section of vmd file to split up sections
    partfile.write(str(len(particle)-1)+"\n"+"break \n")
    # Create loop for each particle
    for f in range(0,len(particle)-1):
        if particle[f] == None:
            break
        else:
            # Collect each planets energy in list
            energy.append(particle[f].kinetic_energy() + potenergy[f])
            # Write VMD file component
            partpos=Particle3D.__str__(particle[f])
            partfile.write(partpos + "\n")
    # Find total energy of system and apply to list & file
    totalenergy.append(sum(energy))
    efile.write(str(sum(energy))+ "\n")


    ############### APO AND PERIAPSES CODE ##################
    # Initialise apo and periapsis lists of required length
    apo = [None] * (len(particle)-1)
    peri = [None] * (len(particle)-1)
    # Initiate calculation of initial apo and periapses
    apo, peri = apoperi(particle, apo, peri)

    ################ MOON SEPARATION ##################
    # Finding the separation of moon to the earth against
    # time to plot to a graph, allows for better understanding
    # of what is an acceptable dt value
    moonsep_list.append(np.linalg.norm(Particle3D.seperation(particle[3].position,particle[4].position)))


    ################### TIME PERIOD CODE ########################
    # IDEA: Take the initial seperation from the sun of each object.
    # Take dot product with new seperation from the sun and measure the times we have a change of sign.
    # Accuracy increases with higher value of time-end as it is averaged over all orbits.
    # A balance between speed and accuracy is chosen.
    p_const_initial = []
    #Initialise counter
    counter = np.zeros((len(particle)-1))
    for i in range(0,len(particle)-1):
        # Excluding moon, relative to sun
        if i!= 4:
            p_const_initial.append(Particle3D.seperation(particle[0].position,particle[i].position))
        # Moon consideration
        else:
            p_const_initial.append(Particle3D.seperation(particle[3].position,particle[i].position))
    ################# TIME INTEGRATION LOOP BEGINS #################
    # Basic time loop, from 0->time_end
    for i in range(0,time_end,dt):

        # Determine new position of each particle
        for f in range(0,len(particle)-1):
            if particle[f] == None:
                break
            else:
                particle[f].step_pos2nd(dt,totalforce[f])

        ######### APO AND PERIAPSES CODE #########
        apo, peri = apoperi(particle, apo, peri)

        ######### TIME PERIOD CODE ######
        for c in range(1,len(particle)-1):
            # Everything relative to Sun
            if c != 4:
                if dot_prod(p_const_initial[c],Particle3D.seperation(particle[0].position,particle[c].position)) <= 0 and particle[c] != particle[0]:
                    counter[c] = counter[c]+0.5 # Because the sign will change about 2 times per rotation
                    p_const_initial[c] = -p_const_initial[c]
            # Moon relative to Earth
            else:
                if dot_prod(p_const_initial[c],Particle3D.seperation(particle[3].position,particle[c].position)) <= 0 and particle[c] != particle[0]:
                    counter[c] = counter[c]+0.5 # Because the sign will change about 2 times per rotation
                    p_const_initial[c] = -p_const_initial[c]


        ###### NEW FORCE AND ENERGY CALCULATION #####
        totalforcenew, potenergy = baseloop(particle)

        # Sets how often the progress is shown to the user, can be edited to suit the user.
        if i % 1000 == 0:
            print("Step " + str(i) + "/" + str(time_end))

        # Determine velocity of each particle
        for f in range(0,len(particle)-1):
            if particle[f] == None:
                break
            else:
                particle[f].step_velocity(dt,0.5*(totalforce[f]+totalforcenew[f])) #0.5


        # Definining new totalforce
        totalforce = totalforcenew

        # Increase in time
        time += dt

        ##################### WRITING TO FILES #####################
        ## This component identical to earlier section, improvement would be to make this into a function.
        # Write total energy to Files
        # Resetting energy value
        energy = []
        partfile.write(str(len(particle)-1)+"\n"+"break \n")
        for f in range(0,len(particle)-1):
            if particle[f] == None:
                break
            else:
                # Write energy of each planet
                energy.append(particle[f].kinetic_energy() + potenergy[f])
                # Write VMD file component
                partpos=Particle3D.__str__(particle[f])
                partfile.write(partpos + "\n")
        # Find total energy of system and apply to list & file
        totalenergy.append(sum(energy))
        efile.write(str(sum(energy))+ "\n")

        ################ MOON SEPARATION ##################
        # Finding the separation of moon to the earth against
        # time to plot to a graph, allows for better understanding
        # of what is an acceptable dt value
        moonsep_list.append(np.linalg.norm(Particle3D.seperation(particle[3].position,particle[4].position)))

        # Append information to data lists
        time_list.append(time)

    ######################## TIME INTEGRATION LOOP ENDS ############################

    # Post-simulation:
    # Make a time correction on the value of the counter.
    # The angle of the last position of the simulation with initial position of each object is divided by the maximum angle attained by our counter system(pi)
    # Gives a first order correction to our method
    for c in range(1,len(particle)-1):
    #Find the Time period of each object based on Earth years.
        # Ensures no division by 0
        if counter[c] != 0:
            # Ensures out timer is averaged across all periods.
            ## Possible issue! If we record 1.4 of a full period, this 0.4 isn't accounted for.
            ## Solved by making the simulation length approx to 2 periods of pluto, as other planets
            ## will average out more and plutos orbit will be correct within error.
            counter[c] = (time_end/365)/counter[c]
    # Close energy and vmd output files
    efile.close()
    partfile.close()



    #Write apo and periapses to files
    apfile = open("apoperifile.txt", "w")
    #creates loop to make file readable while skipping repetition
    for f in range(1,len(particle)-1):
        apfile.write("=============================== \n")
        apfile.write("Particle " + str(particle[f].label) + " - " + str(plabel[f]) + " \n")
        apfile.write("Apoapsis - " + str(format(apo[f],"E")) + " km" " \n" + "Periapsis - " + str(format(peri[f],"E")) + " km" + " \n")
    apfile.write("=============================== \n")

    #close file
    apfile.close()

    #Write period to files
    periodfile = open("time_periodfile.txt", "w")
    for f in range(0,len(particle)-1):
        periodfile.write("=============================== \n")
        periodfile.write("Particle " + str(particle[f].label) + " - " + str(plabel[f]) + (" \n"))
        periodfile.write("Time Period - " + str(counter[f]) +" "+ "Earth years" + " \n")
        periodfile.write("=============================== \n")
    periodfile.close()

    # Plot Total Energy of the system.
    pyplot.title('Solar System simulation: total energy vs time')
    pyplot.xlabel('Time/Days')
    pyplot.ylabel('Energy/1.3996E-4 Joules')
    pyplot.plot(time_list, totalenergy,'r', label ="Total Energy" )
    pyplot.legend(loc = "best")
    pyplot.show()

    # Plot Moon-Earth separation in the system.
    pyplot.title('Moon-Earth Separation: separation vs time')
    pyplot.xlabel('Time/Days')
    pyplot.ylabel('Separation/km')
    pyplot.plot(time_list, moonsep_list,'r', label ="moon position relative to earth" )
    pyplot.legend(loc = "best")
    pyplot.show()

############# RUN CODE ###############

#Load in variables first from file
G,dt,time_end,vmdfile,positions = fileread(sys.argv)
#Execute main code
main()
