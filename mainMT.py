# Joseph Cooney
import math
from numpy.random import rand, normal
import time
import vpython
import threading

import interactions
from oval import Oval
from cone import Cone
from interactions import potLJ, potES, potEwald, getDist
from math import pi, sqrt
from cartPolar import toPolar, toCart


myCanvas = vpython.canvas(background=vpython.color.gray(0.5))
vpython.scene.delete()
TEMP = 260 # Kelvin(s) - 250-273
BOLTZ = 8.617 * (10 ** -5) # eV/K

axes = [0.5859, 0.2, 1.5139]

# min number for realistic ice nucleus is ~8303 molecules
numMols = int(input("enter the number of molecules: ")) # for some reason?
spacing = 5
max_steps = 1000
molecules = []
solidified = []
allMols = []
nSize = 1000 # 1-5000 so that whole nucleus is ~um (10,000 Ang)
# nucleus charge is in interactions file.

for obj in myCanvas.objects:
    obj.delete()

nucleus = Oval([nSize, nSize, nSize], [0, 0, 0], [0, 0, 0], myCanvas)

# bounding parameters
boundPol = pi # theta [0, pi)
boundAz = pi / 3 # phi [0, 2pi)
boundRad = nSize
boundSpeed = 0.1 * sqrt(3 * TEMP * 8.314 / 18) # effectively gives a 0.1ns time-step. (bigger speed = larger time-step)
# vRMS of an ideal gas (18AMU) a 250K is (root(3*250 * 8.314 / 18)) = 18.6m/s = 186 Ang/ns.


def addMol(boundR):
    r = normal(boundR, sqrt(boundR))
    theta = rand() * boundPol
    phi = rand() * boundAz
    rot = rand(3)
    molecules.append(Cone(axes, toCart([r, theta, phi]), rot, myCanvas))


def findEnergy(molecule):
    otherMols = []
    for k in allMols:
        if k != molecule:
            otherMols.append(k)

    moleculePos = molecule.getPos()
    molecule2Pos = [mol2.getPos() for mol2 in otherMols]
    allDists = [getDist(moleculePos, pos2) for pos2 in molecule2Pos]
    minimumDist = min(allDists)
    energy = potLJ(molecule, nucleus) + potES(molecule, nucleus)  # nucleus is "short-range"
    ewalds = []

    for m in range(len(molecule2Pos)):
        # short-range
        if allDists[m] < 8.5:
            energy += potLJ(molecule, otherMols[m]) + potES(molecule, otherMols[m])
        else:
            # PM ewald summation
            ewalds.append(otherMols[m])

    energy += potEwald(molecule, ewalds)

    return energy, minimumDist, molecule2Pos


def simulate(index):
    # choose a molecule
    mol = molecules[index]
    oldDOFs = [mol.getPos(), mol.getRot()]
    newDOFs = rand(6)

    # find current energy
    currE = findEnergy(mol)[0]

    # move and rotate mol here
    newPos = [mol.getPos()[k] + (boundSpeed * (newDOFs[k] - 0.5)) for k in range(3)]
    newRot = [newDOFs[3] * pi, newDOFs[4] * 2 * pi, newDOFs[5] * 2 * pi]
    mol.setPos(newPos[0], newPos[1], newPos[2])
    mol.setRot(newRot[0], newRot[1], newRot[2])

    # find new energy
    newE, minDist, mol2Pos = findEnergy(mol)

    if newE <= currE:
        if minDist < 4:
            solidified.append(mol)
            avgR = sum([toPolar(pos)[0] for pos in mol2Pos]) / len(mol2Pos)
            addMol(avgR)

    elif rand() <= math.pow(math.e, ((currE - newE) / (BOLTZ * TEMP))):
        if minDist < 4:
            solidified.append(mol)
            avgR = sum([toPolar(pos)[0] for pos in mol2Pos]) / len(mol2Pos)
            addMol(avgR)

    else:
        mol.setPos(oldDOFs[0][0], oldDOFs[0][1], oldDOFs[0][2])
        mol.setRot(oldDOFs[1][0], oldDOFs[1][1], oldDOFs[1][2])


def recordData():
    # write data to file
    minDistances = []
    neighborAngles = []
    nucDistances = []

    for k in range(len(allMols)):
        mol = allMols.pop(k)
        positions = [mol2.getPos() for mol2 in allMols]
        pos = mol.getPos()
        tempDists = [getDist(pos, position) for position in positions]
        minDist = min(tempDists)
        neighbor = allMols[tempDists.index(minDist)]

        nucDistances.append(getDist(mol.getPos(), nucleus.getPos()))
        minDistances.append(minDist)  # nearest-neighbor distance
        neighborAngles.append(mol.getAxis().diff_angle(neighbor.getAxis()) / (pi / 180))  # degrees

        allMols.insert(k, mol)

    for l in range(len(minDistances)):
        dataFileNNDist.write(str(round(minDistances[l], 4)) + ",")

    for l in range(len(neighborAngles)):
        dataFileAngle.write(str(round(neighborAngles[l], 4)) + ",")

    for l in range(len(nucDistances)):
        dataFileNucDist.write(str(round(nucDistances[l], 4)) + ",")

    dataFileNNDist.write("\n")
    dataFileNucDist.write("\n")
    dataFileAngle.write("\n")


# prepare files
fileTime = ""
for timePart in time.localtime()[0:6]:
    fileTime += str(timePart) + "-"

fileTime = fileTime[0:-1]

fileNameNNDist = str(numMols) + "mols_" + fileTime + "NNDist.txt"
fileNameNucDist = str(numMols) + "mols_" + fileTime + "NucDist.txt"
fileNameAngle = str(numMols) + "mols_" + fileTime + "Angle.txt"

dataFileNNDist = open(fileNameNNDist, "x")
dataFileNucDist = open(fileNameNucDist, "x")
dataFileAngle = open(fileNameAngle, "x")

dataFileNNDist.write("#Joe Cooney " + fileTime + " ice simulation data for " + str(numMols) + " molecules and "
                     + str(interactions.nCharge) + " nucleus charge\n" + "#neighbor_distance\n")

dataFileNucDist.write("#Joe Cooney " + fileTime + " ice simulation data for " + str(numMols) + " molecules and "
                      + str(interactions.nCharge) + " nucleus charge\n" + "#nucleus_distance\n")

dataFileAngle.write("#Joe Cooney " + fileTime + " ice simulation data for " + str(numMols) + " molecules and "
                    + str(interactions.nCharge) + " nucleus charge\n" + "#neighbor_angle\n")

# summon cones
for i in range(numMols + 1):
    addMol(boundRad)

nThreads = 10

# simulation loop
for i in range(max_steps):
    # startTime = time.time()
    allMols = solidified + molecules
    numMolecules = len(molecules)
    molsPerThread = numMolecules // nThreads
    ind = 0
    allocated = []

    while ind < (numMolecules-nThreads):
        threads = []

        for j in range(nThreads):
            threads.append(threading.Thread(target=simulate, args=(ind + j,)))

        ind += nThreads

        for t in threads:
            t.start()

        for t in threads:
            t.join()

    recordData()

dataFileNNDist.close()
dataFileNucDist.close()
dataFileAngle.close()

# Data analysis: this program saves data as .txt file. Analyze in different program.
# Why? matplotlib is not compatible and I can't downgrade python since this is also incompatible.


if str(input("press enter to end simulation: ")).strip() == "":
    # remove error messages on exit and keep window open
    vpython.sleep(time.time())
