# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 16:03:59 2020

GROUP STUDIES
NA62 EXPERIMENT CODE

@author: Thomas MacFadyen
"""
from ps1 import FourVector, BoostMatrix
from ps2 import ParticleDatabase
import random
import math
import numpy as np
import matplotlib.pyplot as plt

database = ParticleDatabase()

def two_decay(combined_COM_vector, product1, product2):
    """
    Performs a two-body decay in the COM frame and returns the FourVectors of
    the two products.
    By default, these two products are a deuteron and a photon, but they can be
    changed to any particle in the ParticleDatabase from ps2.
    """
    # Finds the invariant mass of the system. In the COM frame, sum(p)=0, so
    # the invariant mass is just the total energy.
    M = combined_COM_vector[0]

    # Sets m1 (the mass of product1). By default, this is the mass of a
    # deuteron, but it can be set to the mass of a photon or a pion by changing
    # the argument of the function.
    m1 = database[product1].mass
    m2 = database[product2].mass

    # Calculates the energies of product1 and product2.
    E1 = (M**2 + m1**2 - m2**2) / (2*M)
    E2 = (M**2 + m2**2 - m1**2) / (2*M)

    # Picks a random theta and phi to determine the direction of product1.
    theta1 = random.random() * 2 * math.pi
    phi1 = random.random() * math.pi

    # Finds each component of momentum for product1 using E^2 = p^2 + m^2 to
    # find the magnitude of p and using spherical polar coordinates to find the
    # components of p.
    p1 = (E1**2 - m1**2) ** 0.5
    p1x = p1 * math.cos(theta1) * math.sin(phi1)
    p1y = p1 * math.sin(theta1) * math.sin(phi1)
    p1z = p1 * math.cos(phi1)

    # Creates and returns FourVectors for the two products of the decay in the
    # COM frame. Note, product2 must have equal and opposite momentum to
    # product1, so explicit calculations of p2x, p2y and p2z are not required.
    product1_COM_vector = FourVector(E1, p1x, p1y, p1z)
    product2_COM_vector = FourVector(E2, -p1x, -p1y, -p1z)
    
    return product1_COM_vector, product2_COM_vector


def three_decay(combined_COM_vector, product1, product2, product3):
    """
    Performs a three-body decay in the COM frame and returns the FourVectors of
    the three products.
    The products can be chosen to be any particle featured in
    ParticleDatabase() from ps2.
    
    Needs to carry out two two-body decays using an intermediate particle p23,
    so the decay goes like p0 ---> p1 + p23, and then p23 ---> p2 + p3. For
    each two-body decay, we only need the masses of the products. 
    """
    # Finds the invariant mass of the system. In the COM frame, sum(p)=0, so
    # the invariant mass is just the total energy.
    M = combined_COM_vector[0]

    products = [product1, product2, product3]
    product_masses = [0, 0, 0]

    # Sets the masses of each product by identifying if it's a deuteron or
    # photon, and if it's neither, finds the particle mass from the database.
    for i in (0,1,2):
        product_masses[i] = database[products[i]].mass
    
    m1 = product_masses[0]
    m2 = product_masses[1]
    m3 = product_masses[2]
    
    # Sets the min and max possible values of m23 and picks a random value for
    # m23 between its min and max values.
    m23_max = M - m1
    m23_min = m2 + m3
    m23 = m23_min + (random.random() * (m23_max - m23_min))
    
    # Finds the momentum of particle 1  and all its components in the COM frame.
    p1 = ((abs((M - m1 - m23) * (M + m1 + m23) * (M + m1 - m23) * (M - m1 + m23)))**0.5) / (2*M)
    theta1 = random.random() * 2 * math.pi
    phi1 = random.random() * math.pi
    p1x = p1 * math.cos(theta1) * math.sin(phi1)
    p1y = p1 * math.sin(theta1) * math.sin(phi1)
    p1z = p1 * math.cos(phi1)
    
    # Calculates the energy of particle 1 and creates a FourVector for it.
    E1 = (p1**2 + m1**2) ** 0.5
    product1_COM_vector = FourVector(E1, p1x, p1y, p1z)
    
    # Calculates the energy of particle 23 considering the magnitude of the
    # momentum of particle 23 will be the same as for particle 1.
    E23 = (p1**2 + m23**2) ** 0.5
    
    # Particle 23 needs to be boosted into its own COM frame (not the COM frame
    # of the whole system) in order to carry out the second decay. In its own
    # COM frame it has no momentum, so M = m23 for the second decay.
    p2 = ((abs((m23 - m2 - m3) * (m23 + m2 + m3) * (m23 + m2 - m3) * (m23 - m2 + m3)))**0.5) / (2*m23)
    theta2 = random.random() * 2 * math.pi
    phi2 = random.random() * math.pi
    p2x = p2 * math.cos(theta2) * math.sin(phi2)
    p2y = p2 * math.sin(theta2) * math.sin(phi2)
    p2z = p2 * math.cos(phi2)
    
    E2 = (p2**2 + m2**2) ** 0.5
    E3 = (p2**2 + m3**2) ** 0.5
    
    product2_vector = FourVector(E2, p2x, p2y, p2z)
    product3_vector = FourVector(E3, -p2x, -p2y, -p2z)
    
    # Boost product 2 and 3 back into the COM frame of the whole system.
    boost_vector = FourVector(E23, p1x, p1y, p1z)
    Boost = BoostMatrix(boost_vector)
    product2_COM_vector = Boost * product2_vector
    product3_COM_vector = Boost * product3_vector
    
    return product1_COM_vector, product2_COM_vector, product3_COM_vector



def decay(incoming_momentum, incoming_particle, product1, product2, product3=None):
    incoming_mass = database[incoming_particle].mass
    incoming_energy = (incoming_momentum**2 + incoming_mass**2)**0.5
    incoming_lab_vector = FourVector(incoming_energy, 0, 0, incoming_momentum)
    
    boost = BoostMatrix(incoming_lab_vector)
    incoming_COM_vector = boost * incoming_lab_vector
    
    return_lab_vector = FourVector(incoming_energy, 0, 0, -1 * incoming_momentum)
    return_lab_boost = BoostMatrix(return_lab_vector)
    
    if product3 == None:
        products = two_decay(incoming_COM_vector, product1, product2)
        product1_lab_vector = return_lab_boost * products[0]
        product2_lab_vector = return_lab_boost * products[1]
        return product1_lab_vector, product2_lab_vector
    
    else:
        products = three_decay(incoming_COM_vector, product1, product2, product3)
        product1_lab_vector = return_lab_boost * products[0]
        product2_lab_vector = return_lab_boost * products[1]
        product3_lab_vector = return_lab_boost * products[2]
        return product1_lab_vector, product2_lab_vector, product3_lab_vector


def simulate(beam_momentum, number_of_particles, beam_particle, product1, product2, product3=None, cutoff=30):
    import time
    start = time.clock()
    particle1 = []
    particle2 = []
    particle3 = []
    for i in range(0, number_of_particles):
        result = decay(beam_momentum, beam_particle, product1, product2, product3)
        particle1.append(result[0])
        particle2.append(result[1])
        try:
            particle3.append(result[2])
        except:
            pass
    
    decayname = r"%s ---> %s + %s" % (str(beam_particle), str(product1), str(product2))
    if product3 != None:
        decayname += " + %s" % product3
        
    p1E = []
    p1x = []
    p1y = []
    p1z = []
    p1p = []
    p1angle = []
    p2E = []
    p2x = []
    p2y = []
    p2z = []
    p2p = []
    p2angle = []
    p3E = []
    p3x = []
    p3y = []
    p3z = []
    p3p = []
    p3angle = []
    
    for i in particle1:
        p1E.append(i[0])
        p1x.append(i[1])
        p1y.append(i[2])
        p1z.append(i[3])
        p = (i[1] ** 2 + i[2] ** 2 + i[3] ** 2) ** 0.5
        pT = (i[1] ** 2 + i[2] ** 2) ** 0.5
        angle = np.arctan(pT/i[3]) * 1000
        p1p.append(p)
        p1angle.append(angle)
    
    for i in particle2:
        p2E.append(i[0])
        p2x.append(i[1])
        p2y.append(i[2])
        p2z.append(i[3])
        p = (i[1] ** 2 + i[2] ** 2 + i[3] ** 2) ** 0.5
        pT = (i[1] ** 2 + i[2] ** 2) ** 0.5
        angle = np.arctan(pT/i[3]) * 1000
        p2p.append(p)
        p2angle.append(angle)
        
    for i in particle3:
        p3E.append(i[0])
        p3x.append(i[1])
        p3y.append(i[2])
        p3z.append(i[3])
        p = (i[1] ** 2 + i[2] ** 2 + i[3] ** 2) ** 0.5
        pT = (i[1] ** 2 + i[2] ** 2) ** 0.5
        angle = np.arctan(pT/i[3]) * 1000
        p3p.append(p)
        p3angle.append(angle)
        
    p1_30 = 0
    p2_30 = 0
    p3_30 = 0
    
    p1or2_30 = 0

    p1or3_30 = 0
    p2or3_30 = 0
    pifany = 0
    
    for p in p1p:
        if p > cutoff:
            p1_30 += 1
    
    for p in p2p:
        if p > cutoff:
            p2_30 += 1
            
    for p in p3p:
        if p > cutoff:
            p3_30 += 1

    
    for i in range(len(p1p)):
        if p1p[i] > cutoff or p2p[i] > cutoff:
            p1or2_30 += 1
        if len(p3p) != 0:
            if p1p[i] > cutoff or p3p[i] > cutoff:
                p1or3_30 += 1
            if p2p[i] > cutoff or p3p[i] > cutoff:
                p2or3_30 += 1
            if p1p[i] > cutoff or p2p[i] > cutoff or p3p[i]>cutoff:
                pifany += 1

    

    
    p1percent = p1_30 / len(p1p)
    p2percent = p2_30 / len(p2p)
    p1or2percent = p1or2_30/len(p1p)
    p1or3percent = p1or3_30/len(p1p)
    p2or3percent = p2or3_30/len(p2p)
    pifanypercent = pifany/len(p1p)
    try:
        p3percent = p3_30 / len(p3p)
    except:
        p3percent = 0
    
    end = time.clock()
    time = end-start
    print("Time taken = %f seconds" % time)
    print(decayname)
    print("The %s has a probability of %f of have a momentum above %f GeV." % (product1, p1percent, cutoff))
    print("The %s has a probability of %f of have a momentum above %f GeV." % (product2, p2percent, cutoff))
    print("The %s has a probability of %f of have a momentum above %f GeV." % (product3, p3percent, cutoff))

    print("The %s or %s has a probability of %f of have a momentum above %f GeV." % (product1, product2, p1or2percent, cutoff))
    print("The %s or %s has a probability of %f of have a momentum above %f GeV." % (product1, product3, p1or3percent, cutoff))
    print("The %s or %s has a probability of %f of have a momentum above %f GeV." % (product2, product3, p2or3percent, cutoff))

    print("The %s or %s or %s has a probability of %f of have a momentum above %f GeV." % (product1, product2, product3, pifanypercent, cutoff))

        
        
    # Max Parameters
    masses = [0,0,0,0]
    masses[0] = database[beam_particle].mass
    masses[1] = database[product1].mass
    masses[2] = database[product2].mass
    try:
        masses[3] = database[product3].mass
    except:
        pass

    lab = FourVector((beam_momentum**2 + masses[0]**2)**0.5,0,0,-1*beam_momentum)
    boost = BoostMatrix(lab)
    
    Emax = []
    Emin = []
    pxmax = []
    pymax = []
    pzmax = []
    pzmin = []
    
    for i in range(1, len(masses)):
        if i == 1:
            a=2
            b=3
        elif i == 2:
            a=1
            b=3
        elif i == 3:
            a=1
            b=2
        restEmax = (masses[0]**2 + masses[i]**2 - (masses[a] + masses[b])**2) / (2 * masses[0])   # Rest Frame
        restpmax = (restEmax**2 - masses[i]**2) ** 0.5   # Rest Frame
        x = boost*FourVector(restEmax, restpmax, 0, 0)
        y = boost*FourVector(restEmax, 0, restpmax, 0)
        z = boost*FourVector(restEmax, 0, 0, restpmax)
        zmin = boost*FourVector(restEmax, 0, 0, -restpmax)
        Emax.append(z[0])
        Emin.append(zmin[0])
        pxmax.append(x[1])
        pymax.append(y[2])
        pzmax.append(z[3])
        pzmin.append(zmin[3])
        
    if product3 != None:
        data = [[p1E, p1x, p1y, p1z, p1p, p1angle], [p2E, p2x, p2y, p2z, p2p, p2angle], [p3E, p3x, p3y, p3z, p3p, p3angle]]
    else:
        data = [[p1E, p1x, p1y, p1z, p1p, p1angle], [p2E, p2x, p2y, p2z, p2p, p2angle]]

    fontsize = 30
    products = [product1, product2, product3]
    
    for particle in data:
        particlenumber = data.index(particle)
        
        plt.figure(figsize=(25,12))
        plt.hist(particle[0], 200)
        plt.axvline(Emax[particlenumber], 0, 999999, color='r')
        plt.axvline(Emin[particlenumber], 0, 999999, color='r')
        plt.title("%s Energy Spectrum" % products[particlenumber].capitalize(), fontsize=fontsize)
        plt.xlabel("Energy (GeV)", fontsize=0.9*fontsize)
        plt.ylabel("Counts", fontsize=0.9*fontsize)
        plt.xticks(fontsize=0.65*fontsize)
        plt.yticks(fontsize=0.65*fontsize)
        plt.plot()

        '''
        plt.figure(figsize=(25,12))
        plt.hist(particle[1], 200)
        plt.axvline(0, 0, 999999, color='k')
        plt.axvline(pxmax[particlenumber], 0, 999999, color='r')
        plt.axvline(-pxmax[particlenumber], 0, 999999, color='r')
        plt.title(r"%s $p_x$ Spectrum" % products[particlenumber].capitalize(), fontsize=fontsize)
        plt.xlabel("x-component of Momentum (GeV/c)", fontsize=0.9*fontsize)
        plt.ylabel("Counts", fontsize=0.9*fontsize)
        plt.xticks(fontsize=0.65*fontsize)
        plt.yticks(fontsize=0.65*fontsize)
        plt.plot()
        
        plt.figure(figsize=(25,12))
        plt.hist(particle[2], 200)        
        plt.axvline(0, 0, 999999, color='k')
        plt.axvline(pymax[particlenumber], 0, 999999, color='r')
        plt.axvline(-pymax[particlenumber], 0, 999999, color='r')
        plt.title(r"%s $p_y$ Spectrum" % products[particlenumber].capitalize(), fontsize=fontsize)
        plt.xlabel("y-component of Momentum (GeV/c)", fontsize=0.9*fontsize)
        plt.ylabel("Counts", fontsize=0.9*fontsize)
        plt.xticks(fontsize=0.65*fontsize)
        plt.yticks(fontsize=0.65*fontsize)
        plt.plot()
        
        plt.figure(figsize=(25,12))
        plt.hist(particle[3], 200)
        plt.axvline(pzmax[particlenumber], 0, 999999, color='r')
        plt.axvline(pzmin[particlenumber], 0, 999999, color='r')
        plt.title(r"%s $p_z$ Spectrum" % products[particlenumber].capitalize(), fontsize=fontsize)
        plt.xlabel("z-component of Momentum (GeV/c)", fontsize=0.9*fontsize)
        plt.ylabel("Counts", fontsize=0.9*fontsize)
        plt.xticks(fontsize=0.65*fontsize)
        plt.yticks(fontsize=0.65*fontsize)
        plt.plot()
        '''
        
        plt.figure(figsize=(25,12))
        plt.hist(particle[4], 200)
        plt.title(r"%s Total Momentum Spectrum" % products[particlenumber].capitalize(), fontsize=fontsize)
        plt.xlabel("Total Momentum (GeV/c)", fontsize=0.9*fontsize)
        plt.ylabel("Counts", fontsize=0.9*fontsize)
        plt.xticks(fontsize=0.65*fontsize)
        plt.yticks(fontsize=0.65*fontsize)
        plt.plot()
        
        plt.figure(figsize=(25,12))
        plt.hist(particle[5], 200)
        plt.title(r"%s Angle Spectrum" % products[particlenumber].capitalize(), fontsize=fontsize)
        plt.xlabel("Opening Angle (mrad)", fontsize=0.9*fontsize)
        plt.ylabel("Counts", fontsize=0.9*fontsize)
        plt.xticks(fontsize=0.65*fontsize)
        plt.yticks(fontsize=0.65*fontsize)
        plt.plot()

if __name__ == '__main__':
    simulate(75, 1000000, 'K+', 'mu+', 'nu_mu', product3=None, cutoff=30)
    simulate(75, 1000000, 'K+', 'pi+', 'pi0', product3=None, cutoff=30)
    simulate(75, 1000000, 'K+', 'pi+', 'pi+', product3='pi-', cutoff=30)
    simulate(75, 1000000, 'K+', 'pi0', 'e+', product3='nu_e', cutoff=30)
    simulate(75, 1000000, 'K+', 'pi0', 'mu+', product3='nu_mu', cutoff=30)
    simulate(75, 1000000, 'K+', 'pi+', 'pi0', product3='pi0', cutoff=30)
    simulate(75, 1000000, 'K+', 'pi+', 'nu_e', product3='nu_ebar', cutoff=30)

'''
if __name__ != '__main__':
    simulate(75, 10, 'K+', 'pi+', 'pi+', product3='pi-', cutoff=30)
'''