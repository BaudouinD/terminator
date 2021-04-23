import numpy as np
import binary_search as bs
import matplotlib.pyplot as plt


class LinearSpline:
    def __init__(self):
        self.tX = []
        self.dirty = False

    def add_entry(self, t, x):
        self.tX.append([t, x])
        self.dirty = True

    def interpolate(self, t):
        if(self.dirty):
            self.Sort()
            self.dirty = False
        y = 0
        abs = [x[0] for x in self.tX]

        if (abs[0] > t):
            return self.tX[0][1]

        elif (abs[len(abs) - 1] < t):
            return self.tX[len(abs) - 1][1]

        else:
            def key(x): return x[0]
            i = bs.search(self.tX, [t, None], key=key)
            if(i == 0):
                t1 = self.tX[1]
            else:
                t1 = self.tX[i-1]
            t2 = self.tX[i]

            return (t2[1] - t1[1]) / (t2[0] - t1[0]) * (t - t1[0]) + t1[1]

        return

    def Sort(self):
        self.tX.sort(key=lambda x: x[0])


class LinearSpline3D:
    def __init__(self):
        self.splineX = LinearSpline()
        self.splineY = LinearSpline()
        self.splineZ = LinearSpline()
        self.initSpline = False

    def add_entry(self, t, x, y, z):
        self.splineX.add_entry(t, x)
        self.splineY.add_entry(t, y)
        self.splineZ.add_entry(t, z)

    def interpolate(self, t):
        return [self.splineX.interpolate(t), self.splineY.interpolate(t), self.splineZ.interpolate(t)]

    def init_spline(self):
        if(self.init_spline):
            largeur = 1.0
            hauteur = 1.0
            pts_mid = [2.0, 0., 0.]
            A = [pts_mid[0], pts_mid[1] - largeur/2, pts_mid[2]]
            B = [pts_mid[0], pts_mid[1], pts_mid[2] + hauteur]
            C = [pts_mid[0], pts_mid[1] + largeur/2, pts_mid[2]]
            self.add_entry(0., A[0], A[1], A[2])
            self.add_entry(1., B[0], B[1], B[2])
            self.add_entry(2., C[0], C[1], C[2])
            self.add_entry(3., A[0], A[1], A[2])
            self.init_spline = True
    
    def walk_trinalg(self, current_position, dest_position):
        if(self.init_spline):
            print(current_position)
            A = current_position
            B = [current_position[0]+(dest_position[0] - current_position[0] )/2, current_position[1]+(dest_position[1] - current_position[1])/2,  0]
            C = dest_position
            print("A=",A)
            print("B=",B)
            print("C=",C)
            self.add_entry(0., A[0], A[1], A[2])
            self.add_entry(1., B[0], B[1], B[2])
            self.add_entry(2., C[0], C[1], C[2])
            self.add_entry(3., B[0], B[1], A[2])
            #self.add_entry(3., (A[0]-C[0])*1/3, (A[1]-C[1])*1/3, A[2])
            #self.add_entry(4., (A[0]-C[0])*2/3, (A[1]-C[1])*2/3, A[2])
            self.init_spline = True
        


if __name__ == "__main__":
    spline = LinearSpline()
    spline.add_entry(0., 0.)
    spline.add_entry(0.5, 0.2)
    spline.add_entry(1.5, -0.4)
    spline.add_entry(2.3, 0.6)

    xs = np.arange(-0.1, 2.5, 0.1)
    ys = []
    for x in xs:
        ys.append(spline.interpolate(x))

    plt.plot(xs, ys)
    plt.show()

    # xs = [0, 1.2, 2.5, 3.3, 4.2, 5]
    # ys = [0, 6, -3, 4, -2, 0]