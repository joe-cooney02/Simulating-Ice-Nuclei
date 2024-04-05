# Joseph Cooney
from math import sin, cos, atan, sqrt


def toCart(crd: list):
    r = crd[0]
    theta = crd[1]
    phi = crd[2]

    return [r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)]


def toPolar(crd: list):
    x = crd[0]
    y = crd[1]
    z = crd[2]

    return [sqrt((x**2 + y**2 + z**2)), atan(sqrt((x**2 + y**2)) / z), atan(y/x)]
