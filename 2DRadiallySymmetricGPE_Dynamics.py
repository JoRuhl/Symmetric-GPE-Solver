import numpy as np
import scipy.linalg as la
from scipy.sparse import diags
import matplotlib.pyplot as plt

##dynamics 
def dynV(x):
  return 0.5*(omega*x)**2 #this is the new post-quench external potential

#load previously found ground state

groundst = np.loadtxt('GroundState_Beta12800_R2_Wall1_M500.csv', dtype='float', delimiter=' ')

#assign columns in dat file to lists
pts = groundst[0,:]
midpts= groundst[1,:]
soln = groundst[2,:]

#parameters
R = 2 # size of computational domain
M = 500 # integer number of grid points to use for calculation
dr = 2*R/(2*M+1) #size of step in spatial discretization


beta=12800
omega = 93.282029652448571
tau = 1e-4 #discretized time step
period = 2*np.pi/(omega*tau)
time = int(100*period)
time_evol = np.empty([time+1,M+1],dtype=complex)
dynPsi_n = soln[:].astype(complex)
time_evol[0] = dynPsi_n


for t in range(time):

  dynPsi_1 = dynPsi_n*np.exp(-0.5*1j*tau*(dynV(midpts)+beta*dynPsi_n*np.conjugate(dynPsi_n)))
  dynPsi_1 =np.append(dynPsi_1[0],dynPsi_1)

  vec_Psi = np.zeros(len(dynPsi_1)-1,dtype='complex')

  for i in range(M):
    vec_Psi[i] = (1/(4*dr*dr*midpts[i]))*(pts[i+1]*dynPsi_1[i+2]+pts[i]*dynPsi_1[i])-((pts[i+1]+pts[i])/(4*dr*dr*midpts[i])+(1j/tau))*dynPsi_1[i+1]
  vec_Psi = np.append([0],vec_Psi)
  vec_Psi[-1]=0

  dynMaindiag = np.empty(M+2,dtype=complex)
  dynUpperdiag = np.empty(M+1, dtype=complex)
  dynLowerdiag = np.empty(M+1,dtype=complex)

  for i in range(M+1):
    dynMaindiag[i] = (pts[i]+pts[i-1])/(4*dr*dr*midpts[i-1])-1j/tau
    dynUpperdiag[i] = -pts[i]/(4*dr*dr*midpts[i-1])
    dynLowerdiag[i]= -pts[i]/(4*dr*dr*midpts[i])
  dynMaindiag[0],dynMaindiag[-1]=1,1
  dynUpperdiag[0]=-1
  dynLowerdiag[-1]=0

  dynMatrix = diags([dynMaindiag,dynUpperdiag,dynLowerdiag],[0,1,-1]).toarray()

  dynPsi_2 = la.solve(dynMatrix,vec_Psi)

  for k in range(len(dynPsi_2)-1):
    dynPsi_n[k] = dynPsi_2[k+1]*np.exp(-0.5*1j*tau*(dynV(midpts[k])+beta*dynPsi_2[k+1]*np.conjugate(dynPsi_2[k+1])))
  time_evol[t+1]= dynPsi_n[:]


np.savetxt('Time_Evol_beta12800_tau1Eminus4_M500.csv',time_evol)

print('complete')
