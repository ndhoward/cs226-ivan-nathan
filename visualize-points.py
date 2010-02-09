import sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math

input = open(sys.argv[1])

# Read in the data
frames = []
indexInFrame = -1
for line in input:
    vals = line.split(',')
    if len(vals) == 3:
        frame = np.empty((3,int(vals[2])))
        T = int(vals[1])
        while len(frames) < (T-1):
            frames.append(None)
        frames.append(frame)
        indexInFrame = 0
    else:
        x = float(vals[0])
        y = float(vals[1])
        z = float(vals[2])
        t = int(vals[3])
        if math.isnan(x) or math.isnan(y) or math.isnan(z):
            print('arg t:' + str(t) + ' index:' + str(indexInFrame))
        frames[t][:,indexInFrame] = [x,y,z]
        indexInFrame = indexInFrame + 1
        
# Visualize the data
#plt.axes(aspect='equal')
fig = plt.figure()
ax = Axes3D(fig)
ax.set_axis_off()
#ax.set_aspect('equal')
#ax.set_autoscale_on(False)
#ax.set_autoscalex_on(False)
#ax.set_autoscaley_on(False)

indx = 0
while frames[indx] == None:
    indx = indx+1

xs = frames[indx][0,:]
ys = frames[indx][1,:]
zs = frames[indx][2,:]
ax.scatter(xs, ys, zs)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show();
