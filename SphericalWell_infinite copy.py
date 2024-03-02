# ------------------------------
# Quantum Potential Well
# ATOMIC UNITS
#
# V=0 for -10<x<10 and psi=0 outside
# ------------------------------
# Finite differences method as developed by Truhlar JCP 10 (1972) 123-132
#
# code by Jordi Faraudo
#
#
import numpy as np
import matplotlib.pyplot as plt

l= 2

#Potential as a function of position, its an spherical potential so follows the next equation
def getV(x):
    potvalue = l*(l+1)/(2*x*x)
    return potvalue

#Discretized Schrodinger equation in n points (FROM 0 to n-1)
def Eq(n,h,x):
    F = np.zeros([n,n])
    for i in range(0,n):
        F[i,i] = -2*((h**2)*getV(x[i]) + 1)
        if i > 0:
           F[i,i-1] = 1
           if i < n-1:
              F[i,i+1] = 1
    return F

#-------------------------
# Main program
#-------------------------
# Interval for calculating the wave function we have a spherical potential with radius L
L = 1
xlower = 0.000001
xupper = 1

#Discretization options
h = 0.001  #discretization in space

#Create coordinates at which the solution will be calculated
x = np.arange(xlower,xupper+h,h)
#grid size (how many discrete points to use in the range from xlower to xupper
npoints=len(x)

print("Using",npoints, "grid points.")

#Calculation of discrete form of Schrodinger Equation
print("Calculating matrix...")
F=Eq(npoints,h,x)

#diagonalize the matrix F
print("Diagonalizing...")
eigenValues, eigenVectors = np.linalg.eig(F)

#Order results by eigenvalue
# w ordered eigenvalues and vs ordered eigenvectors
idx = eigenValues.argsort()[::-1]
w = eigenValues[idx]
vs = eigenVectors[:,idx]

#Energy Level
E = - w/(2.0*h**2)

print("Enter the level of the quantum state:")
o = int(input())

#Print Energy Values
print("RESULTS:")
for k in range(0,o):
  E_exact=(float(k+1)*(np.pi))**2.0/(2.0*L*L)
  print("n=",k,", E(numeric)=%.4f" %E[k],', E(exact)='+'{:.4f}'.format(E_exact))

#Init Wavefunction (empty list with npoints elements)
psi = [None]*npoints

#Calculation of normalised Wave Functions
for k in range(0,len(w)):
	psi[k] = vs[:,k]
	integral = h*np.dot(psi[k],psi[k])
	psi[k] = psi[k]/integral**0.5

#Plot Wave functions
for v in range(0,o):
    plt.plot(x,psi[v],label=r'$\psi_v(x)$, n = ' + str(v))
    plt.title('Spherical potential well Energy')
    plt.legend()
    plt.xlabel(r'$x$ (a.u.)')
    plt.ylabel(r'$\psi(x)$')

plt.show()

print("Bye")
