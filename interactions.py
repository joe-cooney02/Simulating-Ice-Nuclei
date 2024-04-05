# Joseph Cooney
from oval import Oval
from cone import Cone
from math import pow, sqrt, pi
from numpy import fft
from main import nucCharge

# energy scale
E_SCALE = 1.61 # maybe?? This is unclear from the literature!
e0 = 0.005526349406 # e^2 / ev Angstrom
nCharge = nucCharge


# lennard-jones potential between 2 shapes
def potLJ(mol1: Oval or Cone, mol2: Oval or Cone):

    # for oval, this is at center.
    pos1 = mol1.getPos()
    pos2 = mol2.getPos()
    cone1, cone2 = isinstance(mol1, Cone), isinstance(mol2, Cone)

    # for cone, shift to tip
    if cone1:
        p1Shift = mol1.getSites()[0]
        pos1[0] += p1Shift[0]
        pos1[1] += p1Shift[1]
        pos1[2] += p1Shift[2]

    if cone2:
        p2Shift = mol2.getSites()[0]
        pos2[0] += p2Shift[0]
        pos2[1] += p2Shift[1]
        pos2[2] += p2Shift[2]

    # distance between L-J centers - see Abascal, Vega (2005)
    dist = getDist(pos1, pos2)

    if cone1 and cone2:
        # long-range interactions are truncated in main.py
        # Abascal, Vega (2005): 3.1589 angstroms
        size = 3.1589
    else:
        # slightly smaller than it should be by ~0.065 Ang - compensate. (nucleus interactions)
        size = mol1.getAxes()[2] + mol2.getAxes()[2] + 0.065

    return __LJEqn(size, dist)


def __LJEqn(size, dist):
    term = pow(size/dist, 6)
    return 4 * E_SCALE * (pow(term, 2) - term)


def potES(mol1: Oval or Cone, mol2: Oval or Cone):
    mol1Pos, mol2Pos = mol1.getPos(), mol2.getPos()

    pos1 = [mol1Pos]
    pos2 = [mol2Pos]
    charges1, charges2 = [nCharge], [nCharge] # charge of nucleus in e's
    cone1, cone2 = isinstance(mol1, Cone), isinstance(mol2, Cone)

    if cone1:
        p1Shift = mol1.getSites()[1:]
        pos1 = [[mol1Pos[i] + p1Shift[j][i] for i in range(3)] for j in range(3)]  # H, M, H
        charges1 = [0.5564, -1.1128, 0.5564]

    if cone2:
        p2Shift = mol2.getSites()[1:]
        pos2 = [[mol2Pos[i] + p2Shift[j][i] for i in range(3)] for j in range(3)]  # H, M, H
        charges2 = [0.5564, -1.1128, 0.5564]

    uES = 0

    for i in range(len(pos1)):
        for j in range(len(pos2)):
            uES += (charges1[i] * charges2[j]) / (getDist(pos1[i], pos2[j]))

    return (1 / (4 * pi * e0)) * uES


def potEwald(UC: Cone, mols: list[Cone]):
    molCharges = [0.5564, -1.1128, 0.5564]
    charges = [[], [], [], [], []] # charge, distance(s) - cone.pos, then H, M, H
    unitShift = UC.getSites()[1:]
    uPos = [UC.getPos()]
    uPos += [[uPos[0][i] + unitShift[j][i] for i in range(3)] for j in range(3)]

    for mol in mols:
        molPos = mol.getPos()
        shift = mol.getSites()[1:]
        pos = [[molPos[i] + shift[j][i] for i in range(3)] for j in range(3)]

        charges[0].append(0.5564)
        for i in range(1, 5):
            charges[i].append(getDist(pos[0], uPos[i-1]))

        charges[0].append(-1.1128)
        for i in range(1, 5):
            charges[i].append(getDist(pos[1], uPos[i-1]))

        charges[0].append(0.5564)
        for i in range(1, 5):
            charges[i].append(getDist(pos[2], uPos[i-1]))

    pots = [[], []]
    f = 1/(4*pi*e0)
    dists = charges[2:]
    for i in range(len(charges[0])):
        pot = 0
        for c in molCharges:
            for j in range(3):
                pot += (f * c * charges[0][i] / dists[j][i])
        pots[1].append(dists[1])
        pots[0].append(pot)

    chargesFFT = fft.rfft(charges[0:1])
    potsFFT = fft.rfft(pots[0:1])
    return sum([(potsFFT[0][i] * (chargesFFT[0][i]**2)).real for i in range(len(potsFFT[0]))])


def getDist(pos1, pos2):
    return sqrt(pow(pos2[0]-pos1[0], 2) + pow(pos2[1]-pos1[1], 2) + pow(pos2[2]-pos1[2], 2))
