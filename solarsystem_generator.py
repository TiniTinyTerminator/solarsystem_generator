# %%
from collections import defaultdict
import weakref

import plotly.offline as pl
import plotlywidget as plw
import plotly.graph_objs as gro

import pandas as pd

import numpy as np
import numpy.random as rnd

time = 3600 * 24 * 365 * 4
interval = 1800

iterations = int(time/interval)

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter('data.xlsx', engine='xlsxwriter')


def pythagoras3(a, b, c):
    return np.sqrt(a*a+b*b+c*c)


def Gvel(r, m1):
    Gc = 6.67 * pow(10, -11)
    return np.sqrt(Gc * m1/r)


class KeepRefs(object):
    __refs__ = defaultdict(list)

    def __init__(self):
        self.__refs__[self.__class__].append(weakref.ref(self))

    @classmethod
    def get_instances(cls):
        for inst_ref in cls.__refs__[cls]:
            inst = inst_ref()
            if inst is not None:
                yield inst


class spaceObject(KeepRefs):
    Gc = 6.67 * pow(10, -11)

    def __init__(self, name, mass, x=None, y=None, z=None, radius=0.0, inclination=0.0, vStart=0.0, interval=interval):
        if x == None:
            x = rnd.rand()
        if y == None:
            y = rnd.rand()
        if z == None:
            z = rnd.rand()
        super(spaceObject, self).__init__()
        self.active = True
        self.id = rnd.rand()
        self.x = [None for i in range(iterations)]
        self.y = [None for i in range(iterations)]
        self.z = [None for i in range(iterations)]
        self.vx = [None for i in range(iterations)]
        self.vy = [None for i in range(iterations)]
        self.vz = [None for i in range(iterations)]
        self.Fg = [None for i in range(iterations)]

        self.moons = []

        self.Fg[0] = 0.0
        self.interval = interval
        self.name = name
        self.mass = mass
        self.r = radius
        self.x[0] = x
        self.y[0] = y
        self.z[0] = z
        self.vx[0] = 0.0
        self.vy[0] = np.cos(inclination / 180 * np.pi) * vStart
        self.vz[0] = np.sin(inclination / 180 * np.pi) * vStart

    # calculates the Gforce at the given time
    def calcGforce(self, t):
        Fx = 0.0
        Fy = 0.0
        Fz = 0.0

        for j in spaceObject.get_instances():
            if j != self and j.active:
                dist = pythagoras3(
                    j.x[t-1] - self.x[t-1], j.y[t-1] - self.y[t-1], j.z[t-1] - self.z[t-1])

                self.Fg[t] = self.Gc * ((self.mass * j.mass)/dist**2)

                Fx += (self.Fg[t] * (j.x[t-1] - self.x[t-1])) / dist
                Fy += (self.Fg[t] * (j.y[t-1] - self.y[t-1])) / dist
                Fz += (self.Fg[t] * (j.z[t-1] - self.z[t-1])) / dist
        return (Fx, Fy, Fz)

    # calculates the velocity at the given time
    def calcVel(self, t):
        (Fx, Fy, Fz) = self.calcGforce(t)

        Ax = Fx / self.mass
        Ay = Fy / self.mass
        Az = Fz / self.mass

        self.vx[t] = self.vx[t - 1] + Ax * self.interval
        self.vy[t] = self.vy[t - 1] + Ay * self.interval
        self.vz[t] = self.vz[t - 1] + Az * self.interval

        return (self.vx[t], self.vy[t], self.vz[t])

    def calcPos(self, t):
        self.calcVel(t)

        self.x[t] = self.x[t-1] + self.vx[t] * self.interval
        self.y[t] = self.y[t-1] + self.vy[t] * self.interval
        self.z[t] = self.z[t-1] + self.vz[t] * self.interval

        return (self.x[t], self.y[t], self.z[t])

    def addMoon(self, spo, deltaX, deltaY=0, deltaZ=0, deltaV=0, deltaInclination=0, t=0):
        if spo != self and isinstance(spo, spaceObject):
            spo.x[t] = deltaX + self.x[t]
            spo.y[t] = deltaY + self.y[t]
            spo.z[t] = deltaZ + self.z[t]

            v = (deltaV + Gvel(pythagoras3(deltaX, deltaY, deltaZ), self.mass))

            spo.vy[t] = np.cos(deltaInclination / 180 * np.pi) * v + self.vy[t]
            spo.vz[t] = np.sin(deltaInclination / 180 * np.pi) * v + self.vz[t]
            spo.vx[t] = self.vx[t]

    def longestDist(self):
        return np.amax([np.amax(self.x), abs(np.amin(self.x)), np.amax(self.y), abs(np.amin(self.y)), np.amax(self.z), abs(np.amin(self.z))])

    def plot(self, interval=20, Tmin=0, Tmax=None):
        if Tmax != None:
            if Tmax < Tmin:
                Tmax = Tmin + Tmax
            if Tmax >= self.x.__len__():
                Tmax = self.x.__len__() - 1
        else:
            Tmax = self.x.__len__() - 1

        return gro.Scatter3d(x=self.x[Tmin:Tmax:interval], y=self.y[Tmin:Tmax:interval], z=self.z[Tmin:Tmax:interval], name=self.name, mode="lines")

    def plotF(self, interval=3600):
        return gro.Scatter(x=np.arange(1, time, interval), y=self.Fg[1:], name=self.name)

    def plotFG(self, interval=3600):
        arr = []
        for i in range(self.interval, self.Fg.__len__() - 2):
            arr.append(self.Fg[i + 1] - self.Fg[i])

        return gro.Scatter(x=np.arange(1, time, interval), y=arr, name=self.name)

    def render(self, t=None, p=300):
        if t == None or t >= self.x.__len__():
            t = self.x.__len__() - 1

        u = np.linspace(0, np.pi, p)
        v = np.linspace(0, 2 * np.pi, p)

        x = np.outer(np.sin(u), np.sin(v)) * self.r + self.x[t]
        y = np.outer(np.sin(u), np.cos(v)) * self.r + self.y[t]
        z = np.outer(np.cos(u), np.ones_like(v)) * self.r + self.z[t]
        return gro.Surface(x=x, y=y, z=z, name=self.name, showscale=False)

    def toExcel(self):
       # Create a Pandas dataframe from some data.
        df = pd.DataFrame({'x-axis': self.x, 'y-axis': self.y, 'z-axis': self.z,
                           'vel-x': self.vx, 'vel-y': self.vy, 'vel-z': self.vz,
                           'force on object': self.Fg})

        # Convert the dataframe to an XlsxWriter Excel object.
        df.to_excel(writer, sheet_name=self.name)



# generate planets
mercury = spaceObject("mercury", 3.3 * pow(10, 23), radius=2439 * 1000)

venus = spaceObject("Venus", 4.867 * pow(10, 24), radius=6037 * 1000)

earth = spaceObject("earth", 5.972 * pow(10, 24), radius=6371 * 1000)
moon = spaceObject("moon", 7.35 * pow(10, 22), radius=1737 * 1000)
Iss = spaceObject("ISS", 262200, radius=150)

mars = spaceObject("mars", 6.42 * pow(10, 23), radius=6037 * 1000)
phobos = spaceObject("phobos", 	10.6 * pow(10, 15), radius=11000)
deimos = spaceObject("deimos", 	1.4762 * pow(10, 15), radius=6000)

ceres = spaceObject("ceres", 9.393 * pow(10, 20), radius=473 * 1000)

jupiter = spaceObject("jupiter", 1.898 * pow(10, 27), radius=69911 * 1000)
europa = spaceObject("europa", 4800000 * pow(10, 16), radius=3121.6/2 * 1000)
Ganymede = spaceObject("Ganymede", 14819000 *
                       pow(10, 16), radius=5262.4/2 * 1000)
Callisto = spaceObject("Callisto", 10759000 *
                       pow(10, 16), radius=4820.6/2 * 1000)
Io = spaceObject("Io", 10759000 * pow(10, 16), radius=4820.6/2 * 1000)

sun = spaceObject("Sun", 1.989 * pow(10, 30), radius=695508 *
                  1000, vStart=100, inclination=90)

# add set planets positions
sun.addMoon(earth, deltaX=149.6 * pow(10, 9), deltaInclination=7.155)
sun.addMoon(venus, deltaX=108.9 * pow(10, 9), deltaInclination=3.86)
sun.addMoon(mercury, deltaX=57.91 * pow(10, 9), deltaInclination=3.38)
sun.addMoon(mars, deltaX=227.9 * pow(10, 9), deltaInclination=5.65)
sun.addMoon(jupiter, deltaX=778.5 * pow(10, 9), deltaInclination=1.38)
sun.addMoon(ceres, deltaX=414.010 * pow(10, 9), deltaInclination=10)

earth.addMoon(Iss, (408 + 6371) * pow(10, 4))
earth.addMoon(moon, 384400 * pow(10, 3), deltaInclination=5.164)

mars.addMoon(phobos, - 9376 * 1000, deltaInclination=1.093)
mars.addMoon(deimos, 23463.2 * 1000, deltaInclination=0.93)

jupiter.addMoon(europa, 671034 * 1000, deltaInclination=0.471)
jupiter.addMoon(Ganymede, 1070412 * 1000, deltaInclination=0.204)
jupiter.addMoon(Callisto, 1882709 * 1000, deltaInclination=0.205)
jupiter.addMoon(Io, 421700 * 1000, deltaInclination=0.040)

# start calculations
if __name__ == "__main__":

    for j in range(1, iterations):
        for b in spaceObject.get_instances():
            if(j % 100 == 0):
                print(j, '/', iterations, end="\r")
            b.calcGforce(j)
            b.calcVel(j)
            b.calcPos(j)

    objects3d = []
    distances = []
    # start generating 3d-graph
    for i in spaceObject.get_instances():
        objects3d.append(i.plot(interval=300))
        objects3d.append(i.render())
        distances.append(i.longestDist())

    m = np.amax(distances) * 1.002
    pl.plot(dict(data=objects3d,
                 layout=gro.Layout(
                     scene=dict(
                         xaxis=dict(
                             range=[-m, m],
                             title="x-axis [m]"
                         ),

                         yaxis=dict(
                             range=[-m, m],
                             title="y-axis [m]"
                         ),

                         zaxis=dict(
                             range=[-m, m],
                             title="z-axis [m]"
                         )
                     )
                 )
                 ), filename="space.html")
    objects = []
    for i in spaceObject.get_instances():
        objects.append(i.plotFG())
        i.toExcel()

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()
    
    pl.plot(
        dict(
            data=objects,
            layout=gro.Layout(dict(
                scene=dict(
                    yaxis=dict(
                        type='log',
                        title="differential force on object [N]"
                    )
                )
            )
            )
        ), filename="forces.html")
