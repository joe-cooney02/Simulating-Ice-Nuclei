# Joseph Cooney
from vpython import vector, cone, canvas
from math import pi
from cartPolar import toCart


class Cone:
    __axes = vector(0, 0, 0) # length (x, full), height (y, full), width (z, full)
    __pos = vector(0, 0, 0) # x, y, z
    __rot = vector(0, 0, 0) # angles theta, phi, spin orientation of length axis in spherical space.
    __cone = cone(pos=vector(0, 0, 0))
    sites = [vector(0, 0, 0), vector(0, 0, 0), vector(0, 0, 0), vector(0, 0, 0)] # LJ, H, M, H
    __canvas: canvas

    def __init__(self, axes, pos, rot, myCanvas):
        self.__setAxes(axes[0], axes[1], axes[2])
        self.setPos(pos[0], pos[1], pos[2])
        self.setRot(rot[0], rot[1], rot[2])
        axes = toCart([self.__axes.x, self.__rot.x, self.__rot.y])
        self.__cone = cone(pos=self.__pos,
                           axis=vector(axes[0], axes[1], axes[2]),
                           size=vector(self.__axes.x, self.__axes.y, self.__axes.z),
                           canvas=myCanvas,
                           color=vector(0, 0.6, 1))# turns it blue.

        self.__cone.rotate(self.__rot.z)
        self.__canvas = myCanvas
        self.__pos = self.__cone.pos
        self.__setSites()
        self.spinSites(self.__rot.z)

    def __setSites(self):
        self.sites[0] = self.__cone.axis
        tempPos = toCart([self.__axes.z, self.__rot.x, self.__rot.y + (pi/2)])
        self.sites[1] = vector(tempPos[0], tempPos[1], tempPos[2])
        self.sites[2] = self.__cone.axis * 0.8454
        tempPos = toCart([self.__axes.z, self.__rot.x, self.__rot.y - (pi/2)])
        self.sites[3] = vector(tempPos[0], tempPos[1], tempPos[2])

    def spinSites(self, spin):
        self.sites[0] = self.__cone.axis
        self.sites[1].rotate(spin, self.__cone.axis)
        self.sites[2] = self.__cone.axis * 0.8454
        self.sites[3].rotate(spin, self.__cone.axis)

    def __setAxes(self, newX, newY, newZ):
        self.__axes.x = newX
        self.__axes.y = newY
        self.__axes.z = newZ

    def setPos(self, newX, newY, newZ):
        self.__pos.x = newX
        self.__pos.y = newY
        self.__pos.z = newZ
        self.__cone.pos = self.__pos

    def setRot(self, newTheta, newPhi, newSpin):
        self.__rot.x = newTheta
        self.__rot.y = newPhi
        self.__rot.z = newSpin
        axes = toCart([self.__axes.x, self.__rot.x, self.__rot.y])
        self.__cone.axis = vector(axes[0], axes[1], axes[2])
        self.__cone.rotate(self.__rot.z)
        self.spinSites(self.__rot.z)

    def getSites(self):
        s1 = self.sites[0]
        s2 = self.sites[1]
        s3 = self.sites[2]
        s4 = self.sites[3]
        return [[s1.x, s1.y, s1.z], [s2.x, s2.y, s2.z], [s3.x, s3.y, s3.z], [s4.x, s4.y, s4.z]]

    def getAxes(self):
        return [self.__cone.axis.x, self.__cone.axis.y, self.__cone.axis.z]

    def getAxis(self):
        return self.__cone.axis

    def getPos(self):
        return [self.__cone.pos.x, self.__cone.pos.y, self.__cone.pos.z]

    def getRot(self):
        return [self.__rot.x, self.__rot.y, self.__rot.z]

    def getCone(self):
        return self.__cone

    def delete(self):
        self.__cone.visible = False
