import matplotlib.pyplot as plt
import numpy as np
import sys

N = int(sys.argv[1])		# number of bodies
fileName = sys.argv[2]		# name of file containing data

pathToFile = "data/" + fileName + ".txt"
saveFigPositionsPath = "plots/" + fileName + "Pos.png"
saveFigEnergiesPath = "plots/" + fileName + "Erg.png"

# reading body names from first line in the datafile to a list

infile = open(pathToFile, 'r')

earthIndex = 0
listOfBodies = []
for line in infile:
	names = line.split()
	for name in names:
		listOfBodies.append(name)
		if name == "earth":
			earthIndex = names.index(name)		# needed to find max x and y-values later
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
positions = np.zeros((iterations, 2*N))
lineNumber = 0

infile.readline()
for line in infile:
	words = line.split()
	time[lineNumber] = float(words[0])
	kineticEnergy[lineNumber] = float(words[1])
	potentialEnergy[lineNumber] = float(words[2])
	totalEnergy[lineNumber] = float(words[3])
	for i in range(N):
		positions[lineNumber, 2*i] = float(words[2*i + 4])
		positions[lineNumber, 2*i + 1] = float(words[2*i + 5])
	lineNumber += 1

infile.close()

# finding max x and y-values for the earth to check if the orbit is circular
for body in listOfBodies:
	if body == "earth":
		earthXMax = max(positions[:,2*earthIndex])
		earthXMin = min(positions[:,2*earthIndex])
		earthYMax = max(positions[:,2*earthIndex+1])
		earthYMin = min(positions[:,2*earthIndex+1])
		print "earthXMax = %.3f" % earthXMax
		print "earthXMin = %.3f" % earthXMin
		print "earthYMax = %.3f" % earthYMax
		print "earthYMin = %.3f" % earthYMin

# allocating a specific color to each body
listOfColors = []
for body in listOfBodies:
	if body == "sun":
		listOfColors.append('yo')
	elif body == "mercury":
		listOfColors.append('grey')
	elif body == "venus":
		listOfColors.append('lightseagreen')
	elif body == "earth":
		listOfColors.append('blue')
	elif body == "mars":
		listOfColors.append('red')
	elif body == "jupiter":
		listOfColors.append('brown')
	elif body == "saturn":
		listOfColors.append('orangered')
	elif body == "uranus":
		listOfColors.append('green')
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
plt.xlabel('$Distance \ in \ AU$')
plt.ylabel('$Distance \ in \ AU$')
runTime = time[-1]
plt.title('$%.0f \ year \ of \ solar \ system \ evolution$' %runTime)
plt.axis('equal')
plt.legend(listOfBodies)
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
plt.title('$Energy \ of \ the \ Solar \ System$')
plt.legend(['kinetic', 'potential', 'total'])
plt.savefig(saveFigEnergiesPath)

plt.show()

