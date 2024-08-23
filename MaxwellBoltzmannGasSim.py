#GAS SIMULATION WITH MAXWELL-BOLTZMANN STATISTICS

#A simple 2D simulation of colliding circles demonstrating how the Maxwell-Boltzmann distribution of velocities in gasses arises naturally from particle collisions
#After the simulation is complete, a graph is outputted, displaying the experimental distribution of velocities at the end of the sim, comparing them to the theoretical distribution in an ideal gas

import numpy as np
from ipycanvas import Canvas, hold_canvas
from time import sleep
from numba import jit
import matplotlib.pyplot as plt

#Define the canvas size in pixels
w, h = 700, 700  # pixels

canvas = Canvas(width=w, height=h)           #Instantiate an IPyCanvas canvas of height h, width w
display(canvas)                              #Display the canvas
canvas.clear()                               #Clear it
canvas.fill_style = 'white'                  #Set the fill colour to white
canvas.fill_rect(0, 0, width=w, height=h)    #Paint a white rectangle covering the whole canvas
canvas.stroke_rect(0, 0, width=w, height=h)  #Surround the canvas with a rectangle

#Stepper function to update the state of the system (i.e. the position and velocity of the ball)
@jit(nopython=True) #make it speedy
def step(r, v, dt):

#Update the position using the Euler method
    r += v*dt

    #defining the distances between particles
    d2 = (r[jj,0] - r[ii,0])**2 + (r[jj,1] - r[ii,1])**2

    #hitcheck
    hitj, hiti = None, None
    for c in range(npairs):

        if d2[c] < (2*rball)**2:

            hitj, hiti = jj[c], ii[c] #records which two are colliding

            delr = r[hitj,:] - r[hiti,:]#computes positional diference between colliding particles
            delv = v[hitj,:] - v[hiti,:]#computes the velocity difference between colliding particles

            vdotr = delr[0]*delv[0] + delr[1]*delv[1] #computing the dot product of velocity and position

            if vdotr < 0: #checks if particles are moving toward each other

                ssr = delr[0]**2 + delr[1]**2 #the squared distance between the particles

                mj, mi = m[hitj], m[hiti] #finding the masses

                #new bounced velocities
                v[hitj,:] = v[hitj,:] - ((2*mi)/(mi + mj))*((vdotr)/(ssr))*delr
                v[hiti,:] = v[hiti,:] + ((2*mj)/(mi + mj))*((vdotr)/(ssr))*delr





#Computing bouncing off of the walls
    for i in range(N):

        #Left wall
        if (r[i,0] - rball) < 0 and v[i,0] < 0:
            v[i,0] *= -1

        #Right wall
        if (r[i,0] + rball) > w and v[i,0] > 0:
            v[i,0] *= -1

        #Ceiling (which is at height 0 for silly reasons)
        if (r[i,1] - rball) < 0 and v[i,1] < 0:
            v[i,1] *= -1

        #Floor which is at the maximum 'height'
        if (r[i,1] + rball) > h and v[i,1] > 0:
            v[i,1] *= -1

    return r, v

N = 500
rball = 2 #Ball radius

#mass array
m = np.ones(N)


#The state of the particles

#random starting positions
r = np.zeros((N,2))
for i in range(N):
    r[i,:] = [np.random.random()*(w - 2*rball) + rball, np.random.random()*(h - 2*rball) + rball]

#initial velocities of set magnitude and random direction

vinit = 5 #magnitude
v = np.zeros((N,2))

#velocity components are set to be random between -1 and 1, then normalised to a magnitude of 1 and multiplied to the desired magnitude
for i in range(N):
    v[i,:]= [(np.random.random()*2) - 1,(np.random.random()*2) - 1]
    v[i,:]= (v[i,:]/((v[i,0]**2 + v[i,1]**2)**(1/2)))*vinit

jj, ii = np.triu_indices(N, k=1) #defining compressed matrix indices
npairs = len(jj) #number of particle pairs

#energy conservation check
vmag = np.zeros(N, dtype=float)
@jit(nopython=True) #speedy
def meanv(v, vmag):
    speeds = (v[:,0]**2 + v[:,1]**2)**(1/2)
    meanspeed = np.mean(speeds)
    rmsSpeed = (np.mean(speeds**(2)))**(1/2)
    return speeds, meanspeed, rmsSpeed

nsteps = 250
dt = 0.5 #Choose time step
for i in range(nsteps):

    r, v = step(r, v, dt)
    vmag, vav, vrms = meanv(v, vmag)

    with hold_canvas(canvas):
        canvas.clear()
        canvas.fill_style = 'white'
        canvas.fill_rect(0, 0, width=w, height=h)
        canvas.stroke_rect(0, 0, width=w, height=h)
        canvas.fill_style = 'RoyalBlue'
        canvas.fill_circles(r[:,0], r[:,1], rball)
        canvas.fill_style = 'Red'
        canvas.fill_text('i = {:03d}, vav = {:.3f}, vrms = {:04.3f}, vav/vrms = {:04.3f}'.
              format(i, vav, vrms, vav/vrms), 10, h - 10)
    #Sleep command introduces a small delay.
    sleep(0.01)

#plotting the experimental probability density
plt.figure(1, figsize=(12,12))
ax = plt.subplot(111)
pdens = np.histogram(vmag, np.linspace(0, 15, 20, endpoint = True), density=True)
pdens = pdens[:-1] #removing the right hand edge of the last bin
bincentre = np.linspace(0,15,20, endpoint = True)
bincentre = bincentre[:-1] + (bincentre[1] - bincentre[0])/2 #aligning the points with the centre of the bins
ax.plot(bincentre, pdens[0], 'o', 'b') #plotting the points

#plotting the theoretical distribution
v = np.linspace(0, 15, 200, endpoint = True)
mxwll2d = (2/(vrms**2))*v*np.exp(-((v**2)/(vrms**2)))
ax.plot(v, mxwll2d, 'r')

#plotting the experimental speed stats
ax.axvline(vrms, color='k')
ax.axvline(vav, color='Orange')
#ax.axvline(vrms*((np.pi)**(1/2)/2), linestyle='--', color='g')

#labels
ax.legend(['Final speed distribution (simulated)', 'Experimental average velocity', 'Theoretical Maxwell-Boltzmann distribution', 'Experimental RMS velocity'])
ax.set_xlabel('Speed')
ax.set_ylabel('Probability Density (number of particles in speed ranges)')
ax.set_title('Simulated experimental demonstration of the Maxwell-Boltzmann distribution')
