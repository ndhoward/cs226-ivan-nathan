import sys

pointFile = open(sys.argv[1])
output = open('pruned-'+sys.argv[1], 'w')

staticPoints = dict()
numFrames = 0

for line in pointFile:
    vals = line.split(',')
    x = float(vals[0])
    y = float(vals[1])
    z = float(vals[2])
    t = int(vals[3])  
    if t > numFrames:
        numFrames = t

    xIndx = int(x*10)
    zIndx = int(z*10)
    point = frozenset([xIndx,yIndx,zIndx])
    staticPoints[point] = staticPoints.get(point, 0) + 1
pointFile.close()

#Only record points which do not occur at least three times
pointFile = open(sys.argv[1])
oldT = 0
pointList = []
for line in pointFile:
    vals = line.split(',')
    x = float(vals[0])
    y = float(vals[1])
    z = float(vals[2])
    t = int(vals[3])

    if t > oldT:
        output.write('#NewFrame,'+str(oldT)+','+str(len(pointList))+'\n')
        oldT = t
        for line in pointList:
            output.write(line)
        pointList = []
    xIndx = int(x*10)
    yIndx = int(y*10)
    zIndx = int(z*10)
    point = frozenset([xIndx,yIndx,zIndx])
    if staticPoints[point] < 0.1*numFrames:
        pointList.append(line)
output.close()
