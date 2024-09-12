import numpy as np
import scipy.linalg as la
from scipy.sparse import diags
import matplotlib.pyplot as plt

##Find the ground state

def phi(x): # the initial guess for the solution this can be any function
  return np.tanh((1-x)/np.sqrt(np.pi/beta))*np.heaviside(1-x,1)

 

def V(x): #external potential
  return 1e10*np.heaviside(x-wall,1) #the 1 here specifies the value of the Heaviside step function at the wall, otherwise you have a discontinuity

beta = 12800 # dimensionless coupling constant

tau = 1e-4 #discretized time step
tol = 1e-8  #stopping criterion. Difference between updates to solutions

wall = 1 # coordinate of the potential wall if using Heaviside function above, omit if defining different external potential
R = 2 # size of computational domain
M = 500 # integer number of grid points to use for calculation
dr = 2*R/(2*M+1) #size of step in spatial discretization
midpts = np.zeros(M+1) #create empty vectors of the appropriate length to store grid points
pts = np.zeros(M+1)
#define grid points for second order FD
for j in range(M+1): 
  midpts[j] = dr*(j+0.5)

for k in range(M+1): 
  pts[k] = k*dr

# make vector of initial solution guess
soln = np.append([0],[phi(midpts)])
soln[-1] = 0
# make the tridiagonal matrix
maindiag = np.empty(M+2)
upperdiag = np.empty(M+1)
lowerdiag = np.empty(M+1)

for i in range(M+1):
  maindiag[i] = -tau*( (1/(2*dr*dr*midpts[i-1]))*(-pts[i-1]-pts[i])-V(midpts[i-1])-beta*soln[i]*soln[i]-1/tau)
  upperdiag[i] = -tau*pts[i]/(2*dr*dr*midpts[i-1])
  lowerdiag[i] = -tau*pts[i]/(2*dr*dr*midpts[i])

maindiag[0],maindiag[-1] = 1,1
upperdiag[0] = -1
lowerdiag[-1]= 0

matrix = diags([maindiag,upperdiag,lowerdiag],[0,1,-1]).toarray()

#solve first iteration of matrix equation
psi = la.solve(matrix,soln)

# normalize and time advance
diff = 10
iter = 0

while abs(diff) > tol:
  mag = 0
  for p in range(M):
    mag += midpts[p]*psi[p+1]**2

  normsq = 2*np.pi*dr*mag
  diff = np.dot(soln,soln)-np.dot(psi[1:M+2],psi[1:M+2])/normsq
  soln[1:M+1] = psi[1:M+1]/np.sqrt(normsq)
  for i in range(M+1):
    maindiag[i] = -tau*( (1/(2*dr*dr*midpts[i-1]))*(-pts[i-1]-pts[i])-V(midpts[i-1])-beta*soln[i]*soln[i]-1/tau)
    upperdiag[i] = -tau*pts[i]/(2*dr*dr*midpts[i-1])
    lowerdiag[i] = -tau*pts[i]/(2*dr*dr*midpts[i])

  maindiag[0],maindiag[-1] = 1,1
  upperdiag[0] = -1
  lowerdiag[-1]= 0

  matrix = diags([maindiag,upperdiag,lowerdiag],[0,1,-1]).toarray()
  psi = la.solve(matrix,soln)
  iter +=1
  
#save the ground state and computational grid
np.savetxt('GroundState_Beta12800_R2_Wall1_M500.csv',((pts,midpts,soln[1:])))

 #plot the resulting ground state
print('Number of iterations: ', iter)
print('beta:', beta)
plt.plot(pts,soln[1:]*soln[1:])
plt.title('Ground State of 2D Repulsive Condensate in Circular Trap')
plt.xlabel('r')
plt.ylabel('$|\psi_0 (r)|^2$')
plt.show()

