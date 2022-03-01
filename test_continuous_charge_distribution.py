from multipoles import MultipoleExpansion
import numpy as np

# First we set up our grid, a cube of length 10 centered at the origin:

npoints = 101
edge = 10
x, y, z = [np.linspace(-edge/2., edge/2., npoints)]*3
XYZ = np.meshgrid(x, y, z, indexing='ij')


# We model our smeared out charges as gaussian functions:

def gaussian(XYZ, xyz0, sigma):
    g = np.ones_like(XYZ[0])
    for k in range(3):
        g *= np.exp(-(XYZ[k] - xyz0[k])**2 / sigma**2)
    g *= (sigma**2*np.pi)**-1.5
    return g

sigma = 1.5   # the width of our gaussians

# Initialize the charge density rho, which is a 3D numpy array:
rho = gaussian(XYZ, (0, 0, 1), sigma) - gaussian(XYZ, (0, 0, -1), sigma)


# Prepare the charge distribution dict for the MultipoleExpansion object:

charge_dist = {
    'discrete': False,     # we have a continuous charge distribution here
    'rho': rho,
    'xyz': XYZ
}

# The rest is the same as for the discrete case:

l_max = 2   # where to stop the infinite multipole sum; here we expand up to the quadrupole (l=2)

Phi = MultipoleExpansion(charge_dist, l_max)

x, y, z = 30.5, 30.6, 30.7
value = Phi(x, y, z)
print(value)