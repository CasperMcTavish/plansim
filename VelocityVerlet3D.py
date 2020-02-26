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
#For G, in units km^3/(kg*days^2)
#create particle format

##files neeeding loaded in, particle positions, velocities & constants (numstep, dt, G)
##Particle positions taken from 28th February 2020 from NASA

def newton(p1, p2, m1, m2):
    """
    Method to return potential energy and force
    of particle in Newton potential

    :param p1 and p2: Particle3D instances
    :param m1 and m2: masses of the objects under gravitational potential
    :return: potential energy of particle as float & force acting on particle as Numpy array
    """
    sep_scalar_vec = Particle3D.seperation(p1,p2)
    sep_scalar= np.linalg.norm(sep_scalar_vec)
    potential = (G*m1*m2)/sep_scalar
    force = -(G*m1*m2*sep_scalar_vec)/(sep_scalar**3)
    return force, potential



def baseloop(particle):
    """
    Method to process any function across
    all planets and apply results to matrix

    :param particle: Particle list
    :return: Matrix of applied force/energy in suitably sized matrix
    """
    #example of use, totalF = baseloop(particle,force_Newton)
    #Initial conditions

    #Creates w lists, each with h elements.
    #For our use, each list is a particle, each element is a force
    w, h = (len(particle)),(len(particle)-1);
    fmatrix=[[0 for x in range(w)] for y in range(h)]
    pmatrix=[[0 for x in range(w)] for y in range(h)]
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
                #To stop us calculating already existent results
                     if f > g:
                         fmatrix[f][g], pmatrix[f][g] = -fmatrix[g][f], -pmatrix[g][f]
                     else:
                         fmatrix[f][g], pmatrix[f][g] = newton(particle[f].position,particle[g].position, particle[f].mass,particle[g].mass)
                         #print(f, g)
                         #print(element)

    ftotal = [sum(x) for x in fmatrix]
    ptotal = [sum(x) for x in pmatrix]
    return ftotal, ptotal




##########################Begin Main code#####################################
def main():
    # Read name of the two input files from command line
    # Two inputfiles were chosen, one that has information about the particles
    # while the other has information about the interaction and the amount of iterations carried.
    if len(sys.argv) != 3:
        #print("Wrong number of arguments.")
        #print("Usage: " + sys.argv[0] + sys.argv[1] +" Another file")
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
    dt = int(dt)
    time_end = int(time_end)
    #this is only included so I dont have to remove numstep

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
    totalenergy = []

    #Calculate initial total force & total potential energy
    totalforce, potenergy = baseloop(particle)
    #print("Total Forces")
    #print(totalforce)
    #print("Total Energy")
    #print(potenergy)

    #Open FILES
    efile = open("energy.txt","w")
    partfile = open("vmd"+".xyz","w")

    #####################WRITING TO FILES#####################
    #Write total energy to Files
    #Resetting energy value
    #Write total energy to Files
    energy = []
    #Write break section of particle
    partfile.write(str(len(particle)-1)+"\n"+"break \n")
    #Create loop for each particle
    for f in range(0,len(particle)-1):
        if particle[f] == None:
            break
        else:
            #Collect each planets energy in list
            energy.append(particle[f].kinetic_energy() + potenergy[f])
            #Write VMD file component
            partpos=Particle3D.__str__(particle[f])
            partfile.write(partpos + "\n")
    #Find total energy of system and apply to list & file
    totalenergy.append(sum(energy))
    efile.write(str(sum(energy))+ "\n")




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




    # Initialise data lists for plotting later
    #Lists that will contain the time evolution of the position of both particles and their seperation as well as their energy will be taken
    #Since the position of particles is a numpy array, we need to find their norm so we can plot it
    time_list = [time]
    sep_list = [sep]
    energy_list = [totalenergy]

    #################TIME INTEGRATION LOOP#################
    for i in range(0,time_end,dt):
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

        #Update Force, determine new position's potential
        totalforcenew, potenergy = baseloop(particle)

        #Sets how often the progress is shown to the user, can be edited to suit the user.
        if i % 1000 == 0:
            print("Step" + str(i))

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

        #####################WRITING TO FILES#####################
        ##This component identical to earlier section, somehow make this into a function of some sort? Ask on tuesday
        #Write total energy to Files
        #Resetting energy value
        energy = []
        partfile.write(str(len(particle)-1)+"\n"+"break \n")
        for f in range(0,len(particle)-1):
            if particle[f] == None:
                break
            else:
                #Write energy of each planet
                energy.append(particle[f].kinetic_energy() + potenergy[f])
                #Write VMD file component
                partpos=Particle3D.__str__(particle[f])
                partfile.write(partpos + "\n")
        #Find total energy of system and apply to list & file
        totalenergy.append(sum(energy))
        efile.write(str(sum(energy))+ "\n")



        # Append information to data lists
        time_list.append(time)
        #sep_list.append(sep)



    # Post-simulation:
    # Close all output file
    efile.close()
    sepfile.close()
    partfile.close()

    # Plot Total Energy of the system.
    pyplot.title('Velocity Verlet integrator: total energy vs time')
    pyplot.xlabel('Time/Days')
    pyplot.ylabel('Energy/1.3996E-4 Joules')
    pyplot.plot(time_list, totalenergy,'r', label ="Total Energy" )
    pyplot.legend(loc = "best")
    pyplot.show()

# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
