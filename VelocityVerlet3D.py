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



#add the Gravitational constant
#Units used: km, kg and days
G = 4.98217402E-10 
#create particle format

#Creates array divider
#This function is needed so that Forces are added in chunks 
def divide_chunks(l, n):
    """
    Method to divide a list, specifically the force matrix so that, forces are added consistently. 

    :param l: Parameter l is a list.
    :param n: Parameter n is the step size
    :yields: a set fo smaller lists to be added
    """

    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


def force_Newton(p1,p2, m1, m2):
    """
    Method to return the force on a particle
    in a Newtonian potential.

    :param p1 and p2: Particle3D instances
    :param m1: parameter m1 is the mass of the first object
    :param m2: parameter m2 is the mass of the second object
    :return: force acting on particle as Numpy array
    :note this is the Force on p1. To find the force on p2, f1 = -f2 ny Newton's third law
    """
    sep_scalar_vec = Particle3D.seperation(p1,p2)
    sep_scalar= np.linalg.norm(sep_scalar_vec)
    force = -(G*m1*m2*sep_scalar_vec)/(sep_scalar**3)
    return force

#Added m2 section here, this allows us to mirror the arguments of force_Newton, and so can be passed in as a parameter in baseloop
#However m2 is not an active variable
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
                ##Something to improve upon is make this not double count fncs
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
def main():
    # Read name of the two input files from command line
    # Two inputfiles were chosen, one that has information about the particles
    # while the other has information about the interaction and the amount of iterations carried.
    if len(sys.argv) != 3:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + sys.argv[1] +" Another file")
        return
    else:

        inputfile_name1 = sys.argv[1]
        inputfile_name2 = sys.argv[2]

     #Open output files
    positions = open(inputfile_name1, "r")
    int_param = open(inputfile_name2, "r")
    
    # Set up simulation parameters.
    #To portray real physics, the user is asked to input the time the procedure is finished, rather than the amount of steps taken.
    file_handle2 = int_param 
    G,dt,time_end = Particle3D.from_file2(file_handle2)
    G = float(G)
    dt = float(dt)
    time_end = float(time_end)


    
#    # User input is going to be the timestep dt and the amount of times the simulation will be run
#    dt = float(input("The dt is taken to be: "))
#    #To portray real physics, the user is asked to input the time the procedure is finished, rather than the amount of steps taken.
#    t_end = float(input("Please specify the end time as an integer: "))
#    #The number of steps which are needed in the time integration are shown below
    numstep = int(time_end/dt)
    #Initialise time. Note that the natural units given in the question are used. ie timestep is 10.18fs
    time = 0.0



    #Open another text file and extract information about the positiion of the particles and their velocity and creates Particle3D instances stored in list
    file_handle1 = positions
    i=1
#initialisation to be stored in fields.
    particle= []
    particle.append(Particle3D.from_file1(file_handle1))
    #create loop to allow formation of each particle.
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
    totalforce = baseloop(particle,force_Newton)
    print("Total Forces")
    print(totalforce)

  
    #Calculate inital total potential energy
    print("Energies")
    #To create the energy file we need the potential energies of the particle.
    totalenergy = baseloop(particle,pot_energy_Newton) 

    #Write total energy to Files
    energy_list = []
    for f in range(0,len(particle)-1):
        efile = open("energy"+str(f)+".txt","w")
        energy = particle[f].kinetic_energy() + totalenergy[f]
        efile.write(str(energy)+ "\n")
        energy_list.append(energy)

    print("Total Energy")
    #print(totalenergy)
    
    
    #A file or a list for storing the seperation of the two particles is being created and the values are being written into it
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


    # Initialise data lists for plotting later
    #Lists that will contain the time evolution of the position of both particles and their seperation as well as their energy will be taken
    #Since the position of particles is a numpy array, we need to find their norm so we can plot it
    time_list = [time]
    sep_list = [sep]


    #Start the integration loop
    for i in range(numstep):
        #Update the particle position in a loop.
        for f in range(0,len(particle)-1):
            if particle[f] == None:
                break
            else:
                particle[f].step_pos2nd(dt,totalforce[f])

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



        #Write the new positions to VMD file.
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
    # Close all output file
    efile.close()
    sepfile.close()
    partfile.close()

    # Plot Total Energy of the system.
    pyplot.title('Velocity Verlet integrator: total energy vs time')
    pyplot.xlabel('Time/10.18 fs')
    pyplot.ylabel('Energy/eV')
    pyplot.plot(time_list, energy_list,'r', label ="Particle1 Energy" )
    pyplot.legend(loc = "best")
    pyplot.show()
    
# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
