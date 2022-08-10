# abaqus cae noGUI=readODB -- ODB_NAME
# Currently calculates max height and final strain energy

from odbAccess import *
from abaqusConstants import *

from odbMaterial import *
from odbSection import *

import sys

# This is how to print to terminal
# print >> sys.__stdout__, 'Opening ODB'

odbname = sys.argv[-1]

outputDataFile = 'final-' + odbname[:-4] + '-output.txt' # Define output file name

f = open(outputDataFile, "w")


odb = openOdb(odbname)

# Initiate arrays for x, y, z coordinates
xPoss=[]
yPoss=[]
zPoss=[]

lastFrame = odb.steps['Step-Final'].frames[-1] # Read last frame of the simulation
nodalpos = lastFrame.fieldOutputs['COORD'] # Get nodal coordinates
nodalPosvalues = nodalpos.values

for j in range(len(nodalPosvalues)):
	curNodePos = nodalPosvalues[j].data        
	xPoss.append(curNodePos[0])
	yPoss.append(curNodePos[1])                     
	zPoss.append(curNodePos[2])   

ELSE_field  = lastFrame.fieldOutputs['ELSE'].values # Get elemental strain energy
numEl = len(ELSE_field)
total_E = 0
for j in range(numEl):
    total_E = total_E + ELSE_field[j].data

height  = max(zPoss) - min(zPoss)

# Write to file
f.writelines(["%s\n" % height])
f.writelines(["%s\n" % total_E])	

f.close()

odb.close()