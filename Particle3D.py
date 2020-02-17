"""
 CMod Project: Particle3D
 A class to describe 3D particles and their time evolution
 Through the use of the numpy library
 Author: Panayiotis Panayiotou, John Waiton
 StudentID:s1715899

"""
import numpy as np

class Particle3D(object):
    """
    Class to describe 3D particles.

    Properties:
    position(numpy array) - position vector on x,y,z direction
    velocity(numpy array) - velocity vector on x,y,z direction
    mass(integer) - particle mass in amu

    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    * two static methods that reads off values from texts file and implements them as particle instances and characteristics of the interaction
    * static method to calculate the seperation in position of two particle instances
    """
    def __init__(self, label, mass, pos,vel):
        """
        Initialise a Particle3D instance

        :param label: label as string
        :param pos: position as numpy array
        :param vel: velocity as numpy array
        :param mass: mass as integer
        """
        self.position = pos
        self.velocity = vel
        self.mass = mass
        self.label = str(label)


    def __str__(self):
        """
        Define output format.
        The output will be a linewith four parameters: the label, and the x-y-z values of the position vector
        """
        return str(self.label)+" "+ str(self.position[0])+" "+ str(self.position[1]) +" "+str(self.position[2])


    def kinetic_energy(self):
        """
        Calculates the kinetic energy  of 1 of the particles ie the one where our calculations of interactions is based on
        Return kinetic energy as
        1/2*mass*velocity^2
        """
        v = np.linalg.norm(self.velocity)
        v_squared = v**2
        return 0.5*self.mass*v_squared


    def step_velocity(self, dt, force):
        """
        First-order velocity update, as given in the lecture notes
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle as float
        """
        #Since velocity and Force are numpy arrays, we use operation overload for our advantage
        self.velocity += dt*force/self.mass
        #return self.velocity


    def step_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        :note: operation overload can be used but we need to np.add functionality since we were adding a float to an array if we were left with +=
        """
        self.position += (dt*self.velocity)
        #return self.position



    def step_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as float
        """
        self.position += (dt*self.velocity + 0.5*(dt**2)*force/self.mass)
        return self.position


    """
    Static Method 1:
    Seperation distance.
    Computed by finding the difference between the two position vectors

    :param1 p1: The first postion vector
    :param2 p2: The postion vector of the second vector
    """
    @staticmethod
    def seperation(p1,p2):
        displ = p1 - p2
        return displ

    """
    Static Method 2:
    Reading off values from a text file

    :param file_handle: gets values the values of mass position an velocity of the particles interacting
    """
    @staticmethod
    def from_file1(file_handle):
        #read content from file
        #Write the actual code with the linesplit and stuff from last year, will finish at '42' (arbitrary)
        line = file_handle.readline()
        tokens = line.split(" ")
        #creates particle unless list is finished (due to 42 being reached in list)
        if  str(tokens[0]) != "END":
            print(tokens[0] + tokens[1] + tokens[2] + tokens[3])
            pos = np.array([float(tokens[2]),float(tokens[3]),float(tokens[4])])
            vel = np.array([float(tokens[5]),float(tokens[6]),float(tokens[7])])
            #returns Particle
            return Particle3D(float(tokens[0]),float(tokens[1]),pos,vel)
        #if list finished, return empty particle and confirm completion of reading.
        else:
            print("Particle Reading Complete")
            return

    """
    Static Method 2:
    Reading off values from a text file

    :param file_handle: gets values the values of mass position an velocity of the particles interacting
    """
    @staticmethod
    def from_file2(file_handle):
        #read content from file
        #Write the actual code with the linesplit and stuff from last year
        line = file_handle.readline()
        token = line.split(" ")
        De = token[0]
        re = token[1]
        a = token[2]
        return De, re, a
