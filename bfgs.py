import numpy as np
import numpy.linalg as LA
from sympy import symbols, hessian, Function, N, diff

import matplotlib.pyplot as plt
# for line search
# http://www-ai.cs.uni-dortmund.de/LEHRE/VORLESUNGEN/NOPT/SS14/lectures/lecture-07.pdf
# BFGS
# http://terminus.sdsu.edu/SDSU/Math693a/Lectures/18/lecture.pdf

f   = lambda X: 100*(X[1]-X[0]**2)**2 + (1-X[0])**2
g1  = lambda X: -400*X[0]*(-X[0]**2 + X[1]) + 2*X[0] - 2
g2  = lambda X: -200*X[0]**2 + 200*X[1]
h11 = lambda X: 1200*X[0]**2 - 400*X[1] + 2
h12 = lambda X: -400*X[0]
h21 = lambda X: -400*X[0]
h22 = lambda X: 200

def gradient(X):
  g  = np.empty([2,1],dtype=float)
  g[0][0] = g1(X)
  g[1][0] = g2(X)
  return g

def hessian(X):
  H  = np.empty([2,2],dtype=float)
  H[0][0] = h11(X)
  H[0][1] = h12(X)
  H[1][0] = h21(X)
  H[1][1] = h22(X)
  return H

def pk(X):
  H = hessian(X)
  G = gradient(X)
  return np.dot(-1*np.linalg.inv(H),G)

def phi(X,a):
  p = pk(X)
  x_inner = X+a*p
  return f(x_inner).item(0)

def dphi(X, a):
  p = pk(X)
  pT = np.transpose(p)
  x_inner = X+a*p

  G = gradient(x_inner)
  return np.dot(pT, G)

def rosen(x):
    return sum(100*(x[1:]-x[:-1]**2.0)**2.0 +(1-x[:-1])**2.0)

def line_search( X ):
  x_k = X
  c1  = 0.0001
  c2  = 0.9
  a_i = 0
  i   = 1
  phi_0  = phi (x_k, 0)
  dphi_0 = dphi(x_k, 0)
  phi_prev = phi_0
  a_prev   = 0

  #a_i = # SOMETHING HERE for a_1

  while True:
    phi_i = phi(x_k,a_i)              #?? a_i is a_1
    if ( (phi_i > phi_0 + c1*a_i*dphi_0) or ( (phi_i >= phi_prev) and i > 1) ):
      return zoom( a_prev, a_i, c1, c2, X)

    dphi_i = dphi(x_k, a_i) 
    if abs(dphi_i) <= -1*c2*dphi_0:
      return a_i
    if dphi_i >= 0:
      return zoom ( a_i, a_prev, c1, c2, X)
    #a_i = #UPDATE a_i+1 between a_i and a_max
    i+=1
  return

def zoom( a_lo, a_hi, c1, c2, X):
  phi_0   = phi (X, 0)
  dphi_0  = dphi(X, 0)
  a_j     = 0

  while True:
    a_j    = 0.5*( a_lo+a_hi )
    phi_lo = phi(X, a_lo) 
    phi_j  = phi(X, a_j )

    if ( ( phi_j > phi_0 + c1*a_j*dphi_0 ) or (phi_j >= phi_lo) ):
      a_hi = a_j
    else:
      dphi_j = dphi(X, a_j)
      if abs(dphi_j) <= -1*c2*dphi_0 :
        return a_j
      if dphi_j*(a_hi - a_lo) >= 0:
        a_hi = a_lo
      a_lo = a_j
  #return

def bfgs(X, err):
  x_points = []
  y_points = []
  
  I   = np.identity(2)
  x_k = X
  x_points.append(x_k.item(0))
  y_points.append(x_k.item(1))
  H_k = hessian(x_k)
  k = 0
  while abs( LA.norm(gradient(x_k) ) ) > err :
    p_k   = pk(x_k)
    a_k   = line_serach(x_k)
    x_k1  = x_k + a_k*p_k
    s_k   = x_k1 - x_k
    s_kT  = np.transpose(s_k)
    y_k   = gradient(x_k1) - gradient(x_k)
    y_kT  = np.transpose(y_k)
    rho_k = 1./np.dot(y_kT, s_k) 
    H_k1  = np.dot( np.dot ( (I - rho_k * np.dot(s_k, y_kT)) , H_k ) , (I - rho_k * np.dot(y_k,s_kT)) ) 
            + rho_k*np.dot(s_k, s_kT)
    k    += 1

    x_k   = x_k1 
    x_points.append(x_k.item(0))
    y_points.append(x_k.item(1))
    H_k   = H_k1
  return x_k, x_points, y_points

def plot(x_points, y_points):
  xlist = np.linspace(-3,3,100)
  ylist = np.linspace(-3,3,100)

  X,Y = np.meshgrid(xlist,ylist)

  Z = rosen(np.vstack([X.ravel(), Y.ravel()])).reshape((100,100))

  plt.figure() 
  plt.contour(X,Y,Z, np.arange(10)**5)

  # Plotting path of algorithm
  plt.plot(x_points, y_points, '-o')
  plt.show()
  return 

if __name__ == "__main__":
  X             = np.transpose( np.matrix ([[1.2,1.2]]))
  #X             = np.transpose( np.matrix ([[-1.2,1.]]))
  err           = 0.0001

  x_min, x_points, y_points = bfgs(X, err)
  plot(x_points, y_points)
