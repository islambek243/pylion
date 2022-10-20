import pylion as pl
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import numpy as np
from matplotlib.animation import FuncAnimation
# use filename for simulation name
name = Path(__file__).stem

s = pl.Simulation(name)

ions = {'mass': 40, 'charge': -1}
s.append(pl.createioncloud(ions, 1e-3, 1))

trap = {'radius': 3.5e-3, 'length': 2.75e-3, 'kappa': 0.244,'frequency': 4e6, 'voltage': 1000, 'endcapvoltage': 70}
s.append(pl.linearpaultrap(trap))

s.append(pl.langevinbath(1, 1e-5))

s.append(pl.dump('positions.txt', variables=['x', 'y', 'z'], steps=10))
vavg = pl.timeaverage(20, variables=['vx', 'vy', 'vz'])
s.append(pl.dump('secv.txt', vavg, steps=200))

s.append(pl.evolve(1e4))
s.execute()

file1 = open('positions.txt', 'r')
lines = file1.readlines()
coords = []
zdata = np.zeros(1000)
xdata = np.zeros(1000)
ydata = np.zeros(1000)
for i in range(1000):
    coords.append(lines[10*i+9])
    coord = coords[i].rstrip().split()
    xdata[i] = 1e6*float(coord[1])
    ydata[i] = 1e6*float(coord[2])
    zdata[i] = 1e6*float(coord[3])
    #print(str(zdata[i])+'\n')
xav=np.average(xdata)
xdev=np.std(xdata)
yav=np.average(ydata)
ydev=np.std(ydata)
zav=np.average(zdata)
zdev=np.std(zdata)

def func(num, dataSet, lines, redDots):
    # NOTE: there is no .set_data() for 3 dim data...
    line.set_data(dataSet[0:2, :num])
    line.set_3d_properties(dataSet[2, :num])
    redDots.set_data(dataSet[0:2, :num])
    redDots.set_3d_properties(dataSet[2, :num])
    return line


# THE DATA POINTS
t = np.arange(0,1000,1) # This would be the z-axis ('t' means time here)
x = xdata[t]
y = ydata[t]
z = zdata[t]
dataSet = np.array([x, y, z])
numDataPoints = len(t)

# GET SOME MATPLOTLIB OBJECTS
fig = plt.figure()
ax = Axes3D(fig)
redDots = plt.plot(dataSet[0], dataSet[1], dataSet[2], lw=1, c='r', marker='o')[0] # For scatter plot
# NOTE: Can't pass empty arrays into 3d version of plot()
line = plt.plot(dataSet[0], dataSet[1], dataSet[2], lw=1, c='g')[0] # For line plot

# AXES PROPERTIES]
# ax.set_xlim3d([limit0, limit1])
ax.set_xlabel('X(t)')
ax.set_ylabel('Y(t)')
ax.set_zlabel('Z(t)')
ax.set_title('Trajectory of ion for Pauli trap')

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, func, frames=numDataPoints, fargs=(dataSet,lines,redDots), interval=200, blit=False)
line_ani.save(r'Animation.mp4')


plt.show()
