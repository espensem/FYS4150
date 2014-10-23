import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})

N = int(sys.argv[1])		# number of bodies
dt_string = sys.argv[2]
simulationTime_string = sys.argv[3]
fileName = sys.argv[4]		# first part of name of file containing data 
							# ex: SunEarthLocked if name is SunEarthLocked_dt0.001_20yrs
method = sys.argv[5]		# RK4 or Verlet

pathToFile = "data/" + fileName + "_dt" + dt_string + "_" + simulationTime_string + "yrs_" + method + ".txt"
saveFigPositionsPath = "plots/" + fileName + "Pos_dt" + dt_string + "_" + simulationTime_string + "yrs_" + method + ".png"
saveFigEnergiesPath = "plots/" + fileName + "Erg_dt" + dt_string + "_" + simulationTime_string + "yrs_" + method + ".png"
saveFigAngMomPath = "plots/" + fileName + "AngMom_dt" + dt_string + "_" + simulationTime_string + "yrs_" + method + ".png"

dt = float(sys.argv[2])
simulationTime = int(sys.argv[3])
# reading body names from first line in the datafile to a list

infile = open(pathToFile, 'r')

earthIndex = 0
jupiterIndex = 0
listOfBodies = []
for line in infile:
	names = line.split()
	for name in names:
		listOfBodies.append(name)
		if name == "earth":
			earthIndex = names.index(name)		# needed to find max x and y-values later
		if name == "jupiter":
			jupiterindex = names.index(name)	# needed to set plot axis later
	break

infile.close()

# finding the number of iterations

infile = open(pathToFile, 'r')

infile.readline()
iterations = 0
for line in infile:
	iterations += 1

infile.close()

# reading values from file to time vector and position matrix

infile = open(pathToFile, 'r')

time = np.zeros(iterations)
kineticEnergy = np.zeros(iterations)
potentialEnergy = np.zeros(iterations)
totalEnergy = np.zeros(iterations)
angularMomentum = np.zeros(iterations)
positions = np.zeros((iterations, 2*N))
lineNumber = 0

infile.readline()
for line in infile:
	words = line.split()
	time[lineNumber] = float(words[0])
	kineticEnergy[lineNumber] = float(words[1])
	potentialEnergy[lineNumber] = float(words[2])
	totalEnergy[lineNumber] = float(words[3])
	angularMomentum[lineNumber] = float(words[4]) 
	for i in range(N):
		positions[lineNumber, 2*i] = float(words[2*i + 5])
		positions[lineNumber, 2*i + 1] = float(words[2*i + 6])
	lineNumber += 1

infile.close()

# finding max x and y-values for the earth to check if the orbit is circular
for body in listOfBodies:
	if body == "earth":
		earthXMax = max(positions[:,2*earthIndex])
		earthXMin = min(positions[:,2*earthIndex])
		earthYMax = max(positions[:,2*earthIndex+1])
		earthYMin = min(positions[:,2*earthIndex+1])
		print "earthXMax = %.7f" % earthXMax
		print "earthXMin = %.7f" % earthXMin
		print "earthYMax = %.7f" % earthYMax
		print "earthYMin = %.7f" % earthYMin

# allocating a specific color to each body
listOfColors = []
for body in listOfBodies:
	if body == "sun":
		listOfColors.append('bo')
	elif body == "mercury":
		listOfColors.append('grey')
	elif body == "venus":
		listOfColors.append('lightseagreen')
	elif body == "earth":
		listOfColors.append('green')
	elif body == "mars":
		listOfColors.append('brown')
	elif body == "jupiter":
		listOfColors.append('red')
	elif body == "saturn":
		listOfColors.append('orangered')
	elif body == "uranus":
		listOfColors.append('blue')
	elif body == "neptune":
		listOfColors.append('magenta')
	elif body == "pluto":
		listOfColors.append('purple')

# plotting the positions
plt.figure(1)

for i in range(N):
	plt.plot(positions[:, 2*i], positions[:, 2*i + 1], listOfColors[i])
	plt.hold('on')

plt.hold('off')
plt.xlabel('Distance [AU]')
plt.ylabel('Distance [AU]')
#plt.title('$%.0f \ years \ of \ solar \ system \ evolution \ with \ \Delta t=%g$' % (simulationTime, dt))
plt.axis('equal')
plt.legend(listOfBodies)
plt.xlim(-1.6*max(positions[jupiterIndex,:]), 1.6*max(positions[jupiterIndex,:])) 
plt.ylim(-1.2*max(positions[jupiterIndex+1,:]), 1.2*max(positions[jupiterIndex+1,:]))
plt.savefig(saveFigPositionsPath)

# plotting the energy

plt.figure(2)
plt.plot(time, kineticEnergy, 'g')
plt.hold('on')
plt.plot(time, potentialEnergy, 'b')
plt.hold('on')
plt.plot(time, totalEnergy, 'r')
plt.xlabel('$Time \ in \ years$')
plt.ylabel('$Energy \ (dimensionless)$')
#plt.title('$Energy \ of \ the \ Solar \ System \ with \ \Delta t=%g$' % dt)
plt.legend(['kinetic', 'potential', 'total'])
plt.ylim([1.5*min(potentialEnergy), 1.5*max(kineticEnergy)])
plt.savefig(saveFigEnergiesPath)

# plotting the angular momentum

plt.figure(3)
plt.plot(time, angularMomentum)
plt.xlabel('$Time \ in \ years$')
plt.ylabel('$Angular \ momentum \ in \ z-direction \ (dimensionless)$')
#plt.title('$Total \ angular \ momentum \ of \ the \ solar \ system \ with \ \Delta t=%g$' % dt)
plt.savefig(saveFigAngMomPath)

plt.show()

