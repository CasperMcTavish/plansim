"""
Author:Panayiotis Panayiotou and John Waiton
StudentID:S1715899
Updated VelocityVerlet
"""
import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D

#ISSUES - LINE 295 ISH THE FORCES ARE CORRECT BUT ARENT APPLIED CORRECTLY TO THE PARTICLE, VELOCITY IS WRONG

#add the Gravitational constant
#Units used: Astronimical Units(Au), Days, Solar Masses
G = 4.98217402E-10 #Currently in km, kg and days
#2.491627595E-59 for 1
#1.44814E-10 for 0.01
#create particle format

#Creates array divider
def divide_chunks(l, n):

    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]

def force_Newton(p1,p2, m1, m2):
    """
    Method to return the force on a particle
    in a Newtonian potential.

    :param p1 and p2: Particle3D instances
    :param m1: parameter m1 is the mass of the first object
    :param m2: parameter m2 is the mass of the seconf object
    :return: force acting on particle as Numpy array
    :note this is the Force on p1. To find the force on p2, f1 = -f2 ny Newton's third law
    """
    sep_scalar_vec = Particle3D.seperation(p1,p2)
    sep_scalar= np.linalg.norm(sep_scalar_vec)
   # force = (2*De*a*(1-math.exp(-a*(sep_scalar-re)))*(math.exp(-a*(sep_scalar-re))))*(sep_scalar_vec)/sep_scalar
    force = -(G*m1*m2*sep_scalar_vec)/(sep_scalar**3)
    return force

##Added m2 section here, this allows us to mirror the arguments of force_Newton, and so can be passed in as a parameter in baseloop
def pot_energy_Newton(p1, p2, m1, m2):
    """
    Method to return potential energy
    of particle in Newton potential

    :param p1 and p2: Particle3D instances
    :param m1: parameter m1 is the mass of the object under gravitational potential
    :return: potential energy of particle as float
    """
    sep_scalar_vec = Particle3D.seperation(p1,p2)
    sep_scalar= np.linalg.norm(sep_scalar_vec)
    #potential = De*((1-math.exp(math.pow(-a*(sep_scalar-re),2)))-1)
    potential = (G*m1)/sep_scalar
    return potential

def baseloop(particle,fnc):
    """
    Method to process any function across
    all planets and apply results to matrix

    :param particle: Particle list
    :param fnc: Function applied across said Particles, INCLUDING all conditions
    :param matrix: matrix in which results are appended
    :return: Matrix of applied force/energy in suitably size matrix
    """
    #example of use, totalF = baseloop(particle,force_Newton)
    #Initial conditions
    matrix=[]
    #begin basic looping system, allows for looping across all planets
    for f in range(0,len(particle)-1):
        #No particle, break loop. Failsafe
        if particle[f] == None:
            break
        else:
             for g in range(0,len(particle)-1):
                #If we are not testing against ourselves, apply function
                ##Removed the f>g section because it didnt work, and was iterative.
                ##Something to improve upon is make this not double count fncs
                ##Obsolete version of this code available in main()
                 if f != g:
                     if particle[g] == None:
                         break
                     else:
                         element = fnc(particle[f].position,particle[g].position, particle[f].mass,particle[g].mass)
                         print(f, g)
                         print(element)
                         matrix.append(element)

    #Separates array into sections based on size of array
    spl_array = (list(divide_chunks(matrix, len(particle)-2)))
    #Finds the sum of each chunk in x,y,z plane
    total = [sum(l) for l in spl_array]
    return total


##########################Begin Main code#####################################
#I have changed the whole procedure since we only need one file only. since the masses are essentailly in there. IVe Commented most things after checkup remove
def main():
    # Read name of the two input files from command line
    # Two inputfiles were chosen, one that has information about the particle while the other has information about the interaction
    if len(sys.argv) != 2:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] +" Another file")
        return
    else:

        inputfile_name1 = sys.argv[1]
        #inputfile_name2 = sys.argv[2]

     #Open output files
    positions = open(inputfile_name1, "r")
    #potentials = open(inputfile_name2, "r")

    # Set up simulation parameters.
    # User input is going to be the timestep dt and the amount of times the simulation will be run
    dt = float(input("The dt is taken to be: "))
    #To portray real physics, the user is asked to input the time the procedure is finished, rather than the amount of steps taken.
    t_end = float(input("Please specify the end time as an integer: "))
    #The number of steps which are needed in the time integration are shown below
    numstep = int(t_end/dt)
    #Initialise time. Note that the natural units given in the question are used. ie timestep is 10.18fs
    time = 0.0

    #file_handle2 =potentials
    #De, re, g = Particle3D.from_file2(file_handle2)
    #De = float(De)
    #g = float(g)
    #re = float(re)


    #Open another text file and extract information about the positiion of the particle and its velocity and create a Particle3D instance
    #This happens twice
    file_handle1 = positions
    i=1
    #creates false beginning particle to allow initialisation with values 0 0 0 0 0 0 ...
    ##Removed the zeros.
    particle= []
    particle.append(Particle3D.from_file1(file_handle1))
    #create loop to allow formation of each particle with limit of 41 (arbitrary, can be changed easily)
    while True:
        particle.append(Particle3D.from_file1(file_handle1))
        #when particle is empty (list reading finished) end while loop
        if particle[i] == None:
            break
        i+=1
    #show particle number and completion
    print("Particle Processing Complete, particles processed: " + str(len(particle)-1))

    # Get initial force
    #Use the numpy module to create a matrix that will store the values of the force
    totalforce=[]
    totalforcenew = []

    #Calculate initial total force
    print("Forces")
    totalforce = baseloop(particle,force_Newton)


    print("Total Forces")
    print(totalforce)

    #Calculate inital total potential energy
    print("Energies")
    totalenergy = baseloop(particle,pot_energy_Newton)

    #Write total energy to Files
    for f in range(0,len(particle)-1):
        efile = open("energy"+str(f)+".txt","w")
        energy = particle[f].kinetic_energy() + totalenergy[f]
        efile.write(str(energy)+ "\n")

    print("Total Energy")
    print(totalenergy)

    ###########################################################################
    #Total Force Calculations OBSOLETE
    #Adding forces together
    #for f in range(0,len(particle)-1):
        #Take 4 elements from list, E01,E02,etc and adds them together.
    #    g = f*(len(particle)-1)
    #    totalforce.append(np.sum(force_matrix[:(len(particle)-1)],axis = 0))
    ###The following code results in 4 extra files. Resolution same as above
    ###########################################################################

    #To create the energy file we need the potential energies of the particle.  Potiential energy should be positive

    ##First Attempt - TOTAL ENERGY OF PARTICLE 1, 2, 3, 4, 5
    #for f in range(0,len(particle)-1):
    #    if particle[f] == None:
    #        break
    #    else:
    #        efile = open("energy"+str(f)+".txt","w")
    #        energy = particle[f].kinetic_energy() + totalenergy[f]
    #        efile.write(str(energy)+ "\n")

    #A file for dtoring the seperation of the two particles is being created and the values are being written into it
    ##My idea FOR this part of the Code:: for each particle inetraction one file is created ie we should get in our case 16 files. That means it should implement both techniques used above
    ##REMEMBER TO REMOVE - ITERATIVE

    for f in range(0,len(particle)-1):
        if particle[f] == None:
            break
        else:
             for g in range(0,len(particle)-1):
                 #If we are repeating ourselves (f>g), take earlier values and invert them. Otherwise, continue with calculation
                 if particle[g] == None:
                     break
                 elif (f!=g and f<g):
                     sepfile = open("sep"+str(f)+str(g)+".txt","w")
                     sep = np.linalg.norm(Particle3D.seperation(particle[f].position,particle[g].position))
                     sepfile.write(str(sep)+ "\n")

    # Initialise data lists for plotting later
    #Lists that will contain the time evolution of the position of both particles and their seperation as well as their energy will be taken
    #Since the position of particles is a numpy array, we need to find their norm so we can plot it

    #Initialise VMD file creation
    partfile = open("vmd"+".xyz","w")
    partfile.write(str(len(particle)-1)+"\n"+"break \n")
    #Write each initial position of each particle to the file
    for i in range(0,len(particle)-1):
        if particle[i] == None:
            break
        else:
            partpos=Particle3D.__str__(particle[i])
            partfile.write(partpos + "\n")


    time_list = [time]
    sep_list = [sep]
    energy_list = [energy]

    #Write initial positions to file

    #Start the integration loop
    for i in range(numstep):
        #Update the particle position in a loop.
        for f in range(0,len(particle)-1):
            if particle[f] == None:
                break
            else:
                particle[f].step_pos2nd(dt,totalforce[f])
    #wipe the force


        #we dont need any other force
        #Updates Seperation & writes it to file
        for f in range(0,len(particle)-1):
            if particle[f] == None:
                break
            else:
                 for g in range(0,len(particle)-1):
                     #If we are repeating ourselves (f>g), take earlier values and invert them. Otherwise, continue with calculation
                     if particle[g] == None:
                         break
                     elif f != g and f < g:
                         sepfile = open("sep"+str(f)+str(g)+".txt","a")
                         sep = np.linalg.norm(Particle3D.seperation(particle[f].position, particle[g].position))
                         sepfile.write(str(sep) + "\n")

        #Update Force
        print("Forces")
        totalforcenew = baseloop(particle,force_Newton)

        print("Step" + str(i))
        print("Total Forces")
        print(totalforcenew)


        # Determine velocity of each particle
        for f in range(0,len(particle)-1):
            if particle[f] == None:
                break
            else:
                particle[f].step_velocity(dt,0.5*(totalforce[f]+totalforcenew[f])) #0.5


        #Definining new totalforce
        totalforce = totalforcenew

        #Increase in time
        time += dt

        # Determine potential energy
        print("Energies")
        totalenergy = baseloop(particle,pot_energy_Newton)

        #Write total energy to Files
        for f in range(0,len(particle)-1):
            efile = open("energy"+str(f)+".txt","w")
            energy = particle[f].kinetic_energy() + totalenergy[f]
            efile.write(str(energy)+ "\n")

        print("Total Energy")
        print(totalenergy)



        #Write positions to VMD file
        partfile.write(str(len(particle)-1)+"\n"+"break \n")
        #Write each initial position of each particle to the file
        for h in range(0,len(particle)-1):
            if particle[h] == None:
                break
            else:
               partpos=Particle3D.__str__(particle[h])
               partfile.write(partpos + "\n")


        # Append information to data lists
        # Append information to data lists
        time_list.append(time)
        sep_list.append(sep)
        energy_list.append(energy)



    # Post-simulation:
    # Close output file
    efile.close()
    sepfile.close()
    partfile.close()

    # Plot particle trajectory to screen
    #Plotting in one graph, the positions of the particles as well as their seperation
    # The behaviour of the 2 positions should oscillate, as well as their seperation but this is of smaller order that the scales of the graph
    #To avoid this straight line the seperation of the particles has been plotted again
    pyplot.title('Velocity Verlet integrator: position of Particles vs time')
    pyplot.xlabel('Time/10.18 fs')
    pyplot.ylabel('Position/Angstroms')
    pyplot.plot(time_list, pos1_list,'r',label = "Particle1 Position")
    pyplot.plot(time_list, pos2_list,'b',label = "Particle2 Position")

    pyplot.legend(loc = "best")
    pyplot.show()

    # Plot particle energy to screen
    pyplot.title('Velocity Verlet integrator: total energy vs time')
    pyplot.xlabel('Time/10.18 fs')
    pyplot.ylabel('Energy/eV')
    pyplot.plot(time_list, energy_list,'r', label ="Particle1 Energy" )
    pyplot.legend(loc = "best")
    pyplot.show()

    pyplot.title('Velocity Verlet integrator: Seperation of Particles vs time')
    pyplot.xlabel('Time/10.18 fs')
    pyplot.ylabel('Seperation/Angstroms')
    pyplot.plot(time_list, sep_list,'g', label = "Seperation between the Particles")
    pyplot.legend(loc = "best")
    pyplot.show()

# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
