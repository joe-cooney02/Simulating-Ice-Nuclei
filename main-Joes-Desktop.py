# Joseph Cooney
import math
from numpy.random import rand, normal
import time
import vpython
from oval import Oval
from cone import Cone
from interactions import potLJ, potES, getDist
from math import pi, cos, sqrt
from cartPolar import toPolar, toCart


canvas = vpython.canvas(background=vpython.color.gray(0.5))
vpython.scene.delete()
TEMP = 250 # Kelvin(s) - 250-273
BOLTZ = 8.617 * (10 ** -5) # eV/K
STEP_LENGTH = 0.17 # seconds
# this is justification for setting speed to be 0.5A and for preventing it from running too fast.

axes = [0.5859, 0.2, 1.5139]

# min number for realistic ice nucleus is ~8303, nearest cube is 21-> 9261.
numMols = int(input("enter the number of molecules: "))
spacing = 5
max_steps = 1000
molecules = []
solidified = []
allMols = []
nSize = 1000

for obj in canvas.objects:
    obj.delete()

nucleus = Oval([nSize, nSize, nSize], [0, 0, 0], [0, 0, 0], canvas)
solidified.append(nucleus)

# bounding parameters
boundPol = pi / 2 # theta [0, pi)
boundAz = pi / 2 # phi [0, 2pi)
boundRad = nSize


def addMol(boundR):
    r = normal(boundR, sqrt(boundR))
    theta = rand() * boundPol
    phi = rand() * boundAz
    rot = rand(3)
    molecules.append(Cone(axes, toCart([r, theta, phi]), rot, canvas))


# summon cones
for i in range(numMols):
    addMol(boundRad)

# simulation loop
for i in range(max_steps):
    startTime = time.time()
    for j in range(len(molecules)):
        # choose a molecule
        mol = molecules.pop(j)
        oldDOFs = [mol.getPos(), mol.getRot()]
        newDOFs = rand(6)
        allMols = solidified + molecules

        # sum up all energies
        currE = sum([potLJ(mol, mol2) + potES(mol, mol2) for mol2 in allMols])

        # move and rotate mol here
        newPos = [mol.getPos()[i] + (newDOFs[i] - 0.5) for i in range(3)]
        newRot = [newDOFs[3] * pi, newDOFs[4] * 2 * pi, newDOFs[5] * 2 * pi]
        mol.setPos(newPos[0], newPos[1], newPos[2])
        mol.setRot(newRot[0], newRot[1], newRot[2])

        # do PME sum here
        newE = sum([potLJ(mol, mol2) + potES(mol, mol2) for mol2 in allMols])

        molPos = mol.getPos()
        mol2Pos = [mol2.getPos() for mol2 in allMols[1:]]
        minDist = min([getDist(molPos, pos2) for pos2 in mol2Pos])

        if newE <= currE:
            if minDist < 4:
                solidified.append(mol)
                maxR = sum([toPolar(pos)[0] for pos in mol2Pos]) / len(mol2Pos)
                addMol(maxR)
            else:
                molecules.insert(j, mol)

        elif rand() <= math.pow(math.e, ((currE - newE) / (BOLTZ * TEMP))):
            if minDist < 4:
                solidified.append(mol)
                maxR = sum([toPolar(pos)[0] for pos in mol2Pos]) / len(mol2Pos)
                addMol(maxR)
            else:
                molecules.insert(j, mol)

        else:
            mol.setPos(oldDOFs[0][0], oldDOFs[0][1], oldDOFs[0][2])
            mol.setRot(oldDOFs[1][0], oldDOFs[1][1], oldDOFs[1][2])
            molecules.insert(j, mol)

    endTime = time.time()

    if endTime - startTime < STEP_LENGTH:
        time.sleep(STEP_LENGTH - (endTime-startTime))

# do calculations on molecules
if str(input("press enter to calculate density: ")).strip() == "":
    # trace spherical portion containing all molecules, find num molecules (then mass) in the volume
    allMols.pop(0)
    positions = [toPolar(mol.getPos()) for mol in allMols]
    minR = min([pos[0] for pos in positions])
    maxR = max([pos[0] for pos in positions])
    minT = min([pos[1] for pos in positions])
    maxT = max([pos[1] for pos in positions])
    minP = min([pos[2] for pos in positions])
    maxP = max([pos[2] for pos in positions])

    volume = (1/3) * pow(maxR-minR, 3) * (minP-maxP) * (cos(maxT-minT))

    # determine what constitutes a molecule being "in the box?"
    print("density of ice is:", len(allMols) * 18 / volume, "AMU per cubic Angstrom, with", len(allMols), "molecules")

if str(input("press enter to end simulation: ")).strip() == "":
    # remove error messages on exit and keep window open
    vpython.sleep(time.time())
