#COMPRESSIONAL WAVE SIM

#A simple sim of colliding 2D circles with the introduction of a physical perturbation, demonstrating the behaviour of compressional waves (like sound) in fluids
#Once the simulation is complete, a graph is outputted, tracking the populations at the left and right sides of the box, which can be used to demonstrate the speed of 'sound' in the box

import numpy as np
from ipycanvas import Canvas, hold_canvas
from time import sleep
from numba import jit
import matplotlib.pyplot as plt

#Define the canvas size in pixels
w, h = 700, 100  # pixels

leftp = np.zeros(500)
rightp = np.zeros(500)

canvas = Canvas(width=w, height=h)           #Create an IPyCanvas canvas of height h, width w
display(canvas)                              #Display the canvas
canvas.clear()                               #Clear it
canvas.fill_style = 'white'                  #Set the fill colour to white
canvas.fill_rect(0, 0, width=w, height=h)    #Paint a white rectangle covering the whole canvas
canvas.stroke_rect(0, 0, width=w, height=h)  #Surround the canvas with a rectangle

#Stepper function to update the state of the system (i.e. the position and velocity of the ball)
@jit(nopython=True) #makes it speedy
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





#Bounce off walls
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

N = 2000
rball = 2 #Ball radius

#mass array
m = np.ones(N)

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

nsteps = 500
dt = 0.5 #Choose time step
for i in range(nsteps):

    r, v = step(r, v, dt)
    vmag, vav, vrms = meanv(v, vmag)

    for e in range(N):
        if r[e,0] < 21:
            leftp[i] += 1
        if r[e,0] > 679:
            rightp[i] += 1

    #induce perturbation 
    if i == 200:
        for z in range(N):
            if r[z,0] < w/2:
                r[z,0] = w/2

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
    #Sleep command introduces a small delay
    sleep(0.01)






plt.figure(1, figsize=(12,12))
ax = plt.subplot(111)
ax.plot(range(500), leftp, color='blue')
ax.plot(range(500), rightp, color='r')
ax.set_xlabel('time')
ax.set_ylabel('edge populations')
ax.legend(['left side population', 'right side population'])
ax.set_title('Populations of box edges over time')

print('Expected c:')
c = np.sqrt(4/6)*vrms
print(c)


#right hand peaks about 75 steps (150 delta t) after left hand drops
#distance is 350
#speed is about 2.3 p/t, slower than what is expected from the vrms (c = 4.08)
