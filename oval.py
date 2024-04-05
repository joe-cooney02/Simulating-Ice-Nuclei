# Joseph Cooney
from vpython import vpython
from vpython import vector
from math import sin, cos


class Oval:
    # this is used as the nucleus for ice particles
    __axes = vpython.vector(0, 0, 0) # radii in x, y, z
    __pos = vpython.vector(0, 0, 0) # x, y, z
    __rot = vpython.vector(0, 0, 0) # angles theta, phi orientation of major axis in spherical space.
    __oval = vpython.ellipsoid(pos=vector(0, 0, 0))
    canvas: vpython.canvas

    def __init__(self, axes, pos, rot, canvas):
        self.__setAxes(axes[0], axes[1], axes[2])
        self.setPos(pos[0], pos[1], pos[2])
        self.setRot(rot[0], rot[1], rot[2])
        self.__oval = vpython.ellipsoid(pos=self.__pos,
                                        axis=vector((self.__axes.x * sin(self.__rot.x) * cos(self.__rot.y)),
                                                    (self.__axes.x * sin(self.__rot.x) * sin(self.__rot.y)),
                                                    (self.__axes.x * cos(self.__rot.x))),
                                        size=vector(self.__axes.x * 2, self.__axes.y * 2, self.__axes.z * 2),
                                        canvas=canvas)
        self.canvas = canvas
        self.__pos = self.__oval.pos

    def __setAxes(self, newX, newY, newZ):
        self.__axes.x = newX
        self.__axes.y = newY
        self.__axes.z = newZ

    def setPos(self, newX, newY, newZ):
        self.__pos.x = newX
        self.__pos.y = newY
        self.__pos.z = newZ
        self.__oval.pos = self.__pos

    def setRot(self, newTheta, newPhi, newSpin):
        self.__rot.x = newTheta
        self.__rot.y = newPhi
        self.__rot.z = newSpin
        self.__oval.axis = vector((self.__axes.x * sin(self.__rot.x) * cos(self.__rot.y)),
                                  (self.__axes.x * sin(self.__rot.x) * sin(self.__rot.y)),
                                  (self.__axes.x * cos(self.__rot.x)))

    def getAxes(self):
        return [self.__oval.axis.x / 2, self.__oval.axis.y / 2, self.__oval.axis.z / 2]

    def getPos(self):
        return [self.__oval.pos.x, self.__oval.pos.y, self.__oval.pos.z]

    def getRot(self):
        return [self.__rot.x, self.__rot.y, self.__rot.z]

    def getOval(self):
        return self.__oval

    def delete(self):
        self.__oval.visible = False
