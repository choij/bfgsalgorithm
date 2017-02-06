import numpy as np
import numpy.linalg as LA
from sympy import symbols, hessian, Function, N, diff
import matplotlib.pyplot as plt


# http://bicmr.pku.edu.cn/~wenzw/courses/WenyuSun_YaxiangYuan_BB.pdf

f   = lambda X: 100.*(X[1]-X[0]**2)**2 + (1.-X[0])**2
g1  = lambda X: -400.*X[0]*(-X[0]**2 + X[1]) + 2.*X[0] - 2.
g2  = lambda X: -200.*X[0]**2 + 200.*X[1]
h11 = lambda X: 1200.*X[0]**2 - 400.*X[1] + 2.
h12 = lambda X: -400.*X[0]
h21 = lambda X: -400.*X[0]
h22 = lambda X: 200.

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

def bb_stepsize(a_k, g, p_k, x_k):
  s_k  = -(1./a_k)*g
  x_k1 = x_k + s_k
  y_k  = gradient(x_k1) - g

  s_kT = np.transpose(s_k)
  y_kT = np.transpose(y_k)

  xTg  = np.dot(s_kT, y_k)
  gTg  = np.dot(y_kT, y_k)
  print xTg
  print gTg 
  return (xTg)/(gTg)

def grad_desc( X, err ):
  x_points = []
  y_points = []
  x_k = X
  x_points.append(X.item(0))
  y_points.append(X.item(1))
  g   = gradient( x_k )
  p_k = pk( x_k )
  a_k = 1                     # SETTING a_0
  k   = 0

  while abs( LA.norm(a_k*g) ) > err:
      a_k = bb_stepsize(a_k, g, p_k, x_k)

      x_k -= + a_k*g
      x_points.append(X.item(0))
      y_points.append(X.item(1))

      g   = gradient (x_k)
      p_k = pk( x_k )
      k +=1
  print k
  return x_k,x_points,y_points

def pk( X ):
  G = gradient(X)
  Bk = np.identity(2)
  return np.dot(-1.*np.linalg.inv(Bk),G)

if __name__ == "__main__":
  X             = np.transpose( np.matrix ([[1.2,1.2]]))
  #X             = np.transpose( np.matrix ([[-1.2,1.]]))
  err           = 0.00000001

  x_min, x_points, y_points = grad_desc(X, err)
  print x_min
  #plot(x_points, y_points)
