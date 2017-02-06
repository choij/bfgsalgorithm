import numpy as np
import numpy.linalg as LA
from sympy import symbols, Function, N, diff

import matplotlib.pyplot as plt
# for line search
# http://www-ai.cs.uni-dortmund.de/LEHRE/VORLESUNGEN/NOPT/SS14/lectures/lecture-07.pdf
# BFGS
# http://terminus.sdsu.edu/SDSU/Math693a/Lectures/18/lecture.pdf

f   = lambda X: 100.*(X[1]-X[0]**2)**2 + (1.-X[0])**2
g1  = lambda X: -400.*X[0]*(-X[0]**2 + X[1]) + 2.*X[0] - 2.
g2  = lambda X: -200.*X[0]**2 + 200.*X[1]

def gradient(X):
  g  = np.empty([2,1],dtype=float)
  g[0][0] = g1(X)
  g[1][0] = g2(X)
  return g

def pk(X, H):
  G = gradient(X)
  #return np.dot(-1.*np.linalg.inv(H),G)
  return np.dot(-1.*H,G)

def phi(X,a,H):
  p = pk(X, H)
  x_inner = X+a*p
  return f(x_inner).item(0)

def dphi(X, a, H ):
  p = pk(X, H)
  pT = np.transpose(p)
  x_inner = X+a*p

  G = gradient(x_inner)
  return np.dot(pT, G)

def rosen(x):
    return sum(100*(x[1:]-x[:-1]**2.0)**2.0 +(1-x[:-1])**2.0)

def line_search( X, a_max, H ):
  x_k = X
  c1  = 0.0001
  c2  = 0.9
  a_i = 0.
  i   = 1
  phi_0  = phi (x_k, 0., H)
  dphi_0 = dphi(x_k, 0., H)
  phi_prev = phi_0
  a_prev   = a_i

  a_i = (a_prev + a_max)/2.

  while True:
    phi_i = phi(x_k,a_i,H)  
    if ( (phi_i > phi_0 + c1*a_i*dphi_0) or ( (phi_i >= phi_prev) and i > 1) ):
      return zoom( a_prev, a_i, c1, c2, X, H)

    dphi_i = dphi(x_k, a_i,H) 
    if abs(dphi_i) <= -1.*c2*dphi_0:
      return a_i
    if dphi_i >= 0:
      return zoom ( a_i, a_prev, c1, c2, X, H)
    a_i = (a_i + a_max)/2.
    i+=1
  return

def zoom( a_lo, a_hi, c1, c2, X, H):
  phi_0   = phi (X, 0., H)
  dphi_0  = dphi(X, 0., H)
  a_j     = 0.

  while True:
    a_j    = 0.5*( a_lo+a_hi )
    phi_lo = phi(X, a_lo, H) 
    phi_j  = phi(X, a_j, H )

    if ( ( phi_j > phi_0 + c1*a_j*dphi_0 ) or (phi_j >= phi_lo) ):
      a_hi = a_j
    else:
      dphi_j = dphi(X, a_j, H)
      if abs(dphi_j) <= -1.*c2*dphi_0 :
        return a_j
      if dphi_j*(a_hi - a_lo) >= 0.:
        a_hi = a_lo
      a_lo = a_j

def bfgs(X, a_max, err):
  x_points = []
  y_points = []
  
  I   = np.identity(2)
  x_k = X
  x_points.append(x_k.item(0))
  y_points.append(x_k.item(1))

  H_k = I
  k = 0
  while abs( LA.norm(gradient(x_k) ) ) > err :
    p_k   = pk(x_k, H_k)
    a_k   = line_search(x_k, a_max, H_k)
    x_k1  = x_k + a_k*p_k
    s_k   = x_k1 - x_k
    s_kT  = np.transpose(s_k)
    y_k   = gradient(x_k1) - gradient(x_k)
    y_kT  = np.transpose(y_k)
    rho_k = (1./np.dot(y_kT, s_k)).item(0)
    H_k1  = (np.dot( np.dot ( (I - rho_k * np.dot(s_k, y_kT)) , H_k ) , 
            (I - rho_k * np.dot(y_k,s_kT)) ) + rho_k*np.dot(s_k, s_kT))
    k    += 1

    x_k   = x_k1 

    if (k%10 == 0):
      print k, x_k

    x_points.append(x_k.item(0))
    y_points.append(x_k.item(1))
    H_k   = H_k1
  print H_k
  print k
  return x_k, x_points, y_points

def plot(x_points, y_points):
  xlist = np.linspace(-2,2,100)
  ylist = np.linspace(-2,2,100)

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
  err           = 0.001
  a_max         = 1.
  x_min, x_points, y_points = bfgs(X, a_max, err)
  #plot(x_points, y_points)
