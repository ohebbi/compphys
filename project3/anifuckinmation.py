
# Numerical tools
from numpy import *

# Plotting library
from matplotlib.pyplot import *
import mpl_toolkits.mplot3d.axes3d as p3

f = open("values3.txt", "r")
rx = []
ry = []
rz = []
f1 = f.readlines()
for i in f1:
    a = i.split(" ")
    rx.append(float(a[0]))
    ry.append(float(a[1]))
    rz.append(float(a[2]))
f.close()

# Time parameters
T = 1000 # How long to run simulation
dt = 1 # The time step
t = 0
time_steps = int( T / dt ) # Number of time steps

# Put mmatplotlib in interactive mode for animation
ion()

# Setup the figure before starting animation
fig = figure() # Create window
ax = p3.Axes3D(fig)

line, = ax.plot( rx[0], ry[0], rz[0],'bo') # Fetch the line object


# Add other properties to the plot to make it more elegant
fig.suptitle("Jeps") # Title of plot

ax.set_xlim([-2, 2])  # Sets x-axis range
ax.set_ylim([-2, 2])   # Sets y-axis range
ax.set_zlim([-2, 2])
ax.legend(loc='best')   # Adds labels of the lines to the window

draw() # Draws first window
i = 1
# Time loop
while t < T:

    # Update time
    t += dt

    # Plot this new state
    line.set_ydata(ry[i]) # Update the y values of the Psi line
    line.set_xdata(rx[i])
    line.set_zdata(rz[i])
    draw() # Update the plot
    pause(1E-10)
    i += 100

# Turn off interactive mode
ioff()

# Add show so that windows do not automatically close
show()
