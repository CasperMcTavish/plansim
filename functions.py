# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 12:30:10 2020

@author: Panayiotis
"""
import numpy as np
from Particle3D import Particle3D

G = 4.98217402E-10
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


#A function is created to condense the code as this block of code is repeated throughout.
def baseloop(particle):
    """
    Method to process any function across
    all planets and apply results to matrix

    :param particle: Particle list
    :return: Matrix of applied force/energy in suitably sized matrix
    """
    #example of use: totalF = baseloop(particle,force_Newton)
    #Initial conditions

    #Creates w lists, each with h elements. Effectively,creating a list of lists
    #For our use, each list is a particle, each element is a force
    #to initialise our matrices we use list comprehensions.
    w, h = (len(particle)),(len(particle)-1);
    #fmatrix will hold the forces, pmatrix holds the potentials.
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
    #Loop for all planets except sun
    for f in range(1,len(particle)-2):
        #All separations excluding the moon. Special consideration for the moon.
        if f != 4:
            #Calculate separation
            sep = np.linalg.norm(Particle3D.seperation(particle[0].position,particle[f].position))
            #Determine if new apoapsis
            if apo[f] == None:
                apo[f] = sep
            elif (sep > apo[f]):
                apo[f] = sep

            #Determine if new periapsis
            #If no periapsis exists
            if peri[f] == None:
                peri[f] = sep
            elif (sep < peri[f]):
                peri[f] = sep
        #For the moon, relative to the Earth
        else:
            #Calculate separation, assuming Earth is particle 3
            sep = np.linalg.norm(Particle3D.seperation(particle[3].position,particle[f].position))
            #Determine if new apoapsis
            if apo[f] == None:
                apo[f] = sep
            elif (sep > apo[f]):
                apo[f] = sep

            #Determine if new periapsis
            if peri[f] == None:
                peri[f] = sep
            elif (sep < peri[f]):
                peri[f] = sep
    return apo, peri

def dot_prod(a,b):
    dot_prod = np.dot(a,b) 
    value = dot_prod/(np.linalg.norm(a)*np.linalg.norm(b))
    return value
    
