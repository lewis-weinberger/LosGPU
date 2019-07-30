"""
Reads binary output from LosGPU.
Lewis Weinberger 2018
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class losgpu:
 
    def __init__(self, filename):
        """
        Initialises losgpu object by reading the binary output from the
        file at filename. This object has attributes:
            nbins,
            nlos,
            xlos_mid, ylos_mid, zlos_mid,
            xlos, ylos, zlos,
            rho,
            T,
            v,
        which contain the sightline data.
        """
        with open(filename, 'rb') as f:
            # Number of bins and number of sightlines
            self.nbins   = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.nlos    = np.fromfile(f, dtype=np.int32, count=1)[0]
            self.lenlos  = np.fromfile(f, dtype=np.float64, count=1)[0] # cMpc/h

            # Number of line segments in each sightline
            self.nseg  = np.fromfile(f, dtype=np.int32, count=self.nlos) 

            # Coordinates of midpoints for each sightline
            self.xlos_mid  = np.fromfile(f, dtype=np.float64, count=self.nlos) 
            self.ylos_mid  = np.fromfile(f, dtype=np.float64, count=self.nlos)
            self.zlos_mid  = np.fromfile(f, dtype=np.float64, count=self.nlos)

            # Coordinates for each sightline
            self.xlos = np.fromfile(f, dtype=np.float64,
                                     count=self.nbins*self.nlos) # ckpc/h
            self.ylos = np.fromfile(f, dtype=np.float64,
                                     count=self.nbins*self.nlos) # ckpc/h
            self.zlos = np.fromfile(f, dtype=np.float64,
                                     count=self.nbins*self.nlos) # ckpc/h

            # Density, velocity and temperature data
            self.rho = np.fromfile(f, dtype=np.float64,
                                     count=self.nbins*self.nlos) # [g cm^-3]
            self.v   = np.fromfile(f, dtype=np.float64,
                                     count=self.nbins*self.nlos) # [km/s]
            self.T   = np.fromfile(f, dtype=np.float64,
                                     count=self.nbins*self.nlos) # [K]

        print("Sightlines loaded from: ", filename)
        print("NUMLOS = ", self.nlos)
        print("NBINS = ", self.nbins)
        print("Mean number of sightline segments = ", np.mean(self.nseg))

        # Useful quantities
        self.cell_size       = self.lenlos / self.nbins              # cMpc/h
        self.pixel_positions = np.linspace(0, self.lenlos*1e3, 
                                           num=self.nbins)           # ckpc/h

    def plot_los(self, index1, index2, filename):
        """
        Simply plotting functionality for viewing sightline data.
        """
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

        ax1.plot(self.rho[index1*self.nbins:(index1+1)*self.nbins], c='r')
        ax2.plot(self.v[index1*self.nbins:(index1+1)*self.nbins], c='r')
        ax3.plot(self.T[index1*self.nbins:(index1+1)*self.nbins], c='r')
        
        ax1.plot(self.rho[index2*self.nbins:(index2+1)*self.nbins], c='b')
        ax2.plot(self.v[index2*self.nbins:(index2+1)*self.nbins], c='b')
        ax3.plot(self.T[index2*self.nbins:(index2+1)*self.nbins], c='b')

        ax3.set_xlim(0, self.nbins)
        ax3.set_xlabel("Pixel")

        ax1.set_ylabel("density contrast")
        ax2.set_ylabel("v [km s^-1]")
        ax3.set_ylabel("T [K]")

        plt.savefig(filename)

    def plot_los_coords(self, index1, index2, filename, boxlen=160, theta=None, phi=None):
        """
        Simply plotting functionality for viewing sightline position.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        elev, azim = ax.elev, ax.azim

        if phi != None and theta != None:
            ax.view_init(theta, phi)
        elif phi != None and theta == None:
            ax.view_init(elev, phi)
        elif phi == None and theta != None:
            ax.view_init(theta,azim)

        ax.set_xlim(0, boxlen)
        ax.set_ylim(0, boxlen)
        ax.set_zlim(0, boxlen)

        x = self.xlos[index1*self.nbins:(index1+1)*self.nbins]*1e-3
        y = self.ylos[index1*self.nbins:(index1+1)*self.nbins]*1e-3
        z = self.zlos[index1*self.nbins:(index1+1)*self.nbins]*1e-3
        ax.plot(x,y,z,'r.')
        ax.plot([self.xlos_mid[index1]], 
                [self.ylos_mid[index1]], 
                [self.zlos_mid[index1]], 'ro')
        
        x = self.xlos[index2*self.nbins:(index2+1)*self.nbins]*1e-3
        y = self.ylos[index2*self.nbins:(index2+1)*self.nbins]*1e-3
        z = self.zlos[index2*self.nbins:(index2+1)*self.nbins]*1e-3
        ax.plot(x,y,z,'b.')
        ax.plot([self.xlos_mid[index2]], 
                [self.ylos_mid[index2]], 
                [self.zlos_mid[index2]], 'bo')

        ax.set_xlabel("xlos [cMpc/h]")
        ax.set_ylabel("ylos [cMpc/h]")
        ax.set_zlabel("zlos [cMpc/h]")

        plt.savefig(filename)
