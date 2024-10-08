#GRAVITY SIM
#A simple 2D sim of many bodies acting upon each other through gravity

import numpy as np
from ipycanvas import Canvas, hold_canvas
from time import sleep
import matplotlib.pyplot as plt
#Define the window size in pixels
size = w, h = 700, 700 #pixels

canvas = Canvas(width=w, height=h)
display(canvas)
canvas.clear()
canvas.fill_style = 'black'
canvas.fill_rect(0, 0, width=w, height=h)

l, r, b, t = -3E17, 3E17, -3E17, 3E17 #define window bounds in sim units: left, right, bottom, top
rw = r - l                  #window width in metres
rh = t - b                  #window height in metres
pxscl = w/rw                #pixels per m

#Convert from real (simulation) space to window pixel using the scaling value above
def rtopix(x, y):
    pixx = (x - l)*pxscl
    pixy = h - (y - b)*pxscl # y is reversed (the origin is the upper left corner)
    return pixx, pixy

nstars = 100

x = np.random.rand(nstars) - 0.5 #random unscaled + and - positions
x = 2*x #scaled to vary between 1 and -1
x = x*(rw/6) #scales the positions so they span a third of the window in each direction

#same for y
y = np.random.rand(nstars) - 0.5
y = 2*y
y = y*(rw/6)

#defining initial velocities
vx, vy = np.zeros(len(x)), np.zeros(len(y))


#defining star properties
rstar = 1E10
mstar = 2E30

G = 6.67384E-11

zm = 1E6 #zoom factor, originally 1E6

e = 2E16 #softening potential (stops absurd a near 0 distance)


#Update function
@jit(nopython=True) #speedy
def step(x, y, vx, vy, dt):
    """This function advances the animation one time step"""


    #defining the distance between every star
    delx = x[:, np.newaxis] - x #creates an array such that x[a, b] is star a position - star b position, so is the distance b -> a
    dely = y[:, np.newaxis] - y


    dsquared = delx**2 + dely**2 + e**2


    #Define the acceleration
    ax = ((-G*mstar)/((np.sqrt(dsquared))**3))*delx
    ay = ((-G*mstar)/((np.sqrt(dsquared))**3))*dely

    ax = np.sum(ax, axis=1)
    ay = np.sum(ay, axis=1) #collapses the aceleration arrays into 1 dimension, giving the summed accelerations experienced by every star

    #Update the positions and velocities
    vx += ax*dt
    vy += ay*dt
    x  += vx*dt
    y  += vy*dt

    return x, y, vx, vy

dt = 3.1536E7*2E5 #Time step 100,000 hy
nsteps = 1000

for i in range(nsteps):


    x, y, vx, vy = step(x, y, vx, vy, dt)
    pxx, pxy = rtopix(x, y)

    with hold_canvas(canvas):
        canvas.clear()
        canvas.fill_style = 'black'
        canvas.fill_rect(0, 0, width=w, height=h)
        canvas.fill_style = 'LightGray'
        canvas.fill_circles(pxx, pxy, rstar*zm*pxscl)

     #Sleep command introduces a small delay.
    sleep(0.01)
