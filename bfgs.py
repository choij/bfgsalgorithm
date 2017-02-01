import numpy as np
import numpy.linalg as LA
from sympy import symbols, hessian, Function, N, diff

import matplotlib.pyplot as plt

def rosenbrock():
  x, y = symbols('x y')
  f    = symbols('f', cls=Function)
  f    = 100*(y - x**2)**2 + (1 - x)**2
  J11  = diff(f,x)
  J12  = diff(f,y)
  J    = [J11,J12]
  H    = hessian(f, [x,y])
  return f,J,H,x,y

def phi(f,x1,x2,X,a):
  p_i     =
  x_inner = X+a*p_k
  return f.subs([(x1,X[0][0]),(x2,X[1][0])])

def evalu(f,x1,x2,X):
  return f.subs([(x1,X[0][0]),(x2,X[1][0])])

def rosen(x):
    return sum(100*(x[1:]-x[:-1]**2.0)**2.0 +(1-x[:-1])**2.0)

def pk(J, H, x1, x2, X, method):
  g11 = J[0].subs([(x1,X[0]),(x2,X[1])])
  g21 = J[1].subs([(x1,X[0]),(x2,X[1])])

  g  = np.empty([2,1],dtype=float)
  g[0][0] = g11
  g[1][0] = g21

  if (method == 0):
    Bk = np.identity(2)
  else:
    B = H.subs([(x1,X[0]),(x2,X[1])])
    Bk = np.zeros((2,2))
    Bk[0][0] = B[0]
    Bk[0][1] = B[1]
    Bk[1][0] = B[2]
    Bk[1][1] = B[3]
  return np.dot(-1*np.linalg.inv(Bk),g), g

def line_search( f, J, H, x1,x2, X, rho, c, method=0):
  x_k = X
  c1 = 0.0001
  c2 = 0.9
  a_i = 0
  i   = 1
  phi_0  = evalu(phi, 0) #??
  dphi_0 = evalu(dphi, 0)
  phi_prev = phi_0
  a_prev   = a_i

  # set a_1
  while True:
    phi_i = evalu(phi, a_i) #?? for i = 1
    if ( (phi_i > phi_0 + c1*a_i*dphi_0) or ( (phi_i >= phi_prev) and i > 1) ):
      return zoom ( ### a_prev, a_i ... ### )

    dphi_i = evalu(dphi, a_i) #??
    if abs(dphi_i) <= -1*c2*dphi_0:
      return a_i
    if dphi_i >= 0:
      return zoom ( ### a_i, a_prev ... ### )
    a_i = #UPDATE a_i+1 between a_i and a_max
    i+=1
  return

def zoom( a_lo, a_hi, phi, dphi, c1, c2):
  phi_0   = evalu (phi, 0)  #??
  dphi_0  = evalu (dphi, 0) #??
  a_j     = 0

  while True:
    a_j    = 0.5*( a_lo+a_hi )
    phi_lo = evalu( phi, a_lo) #?? 
    phi_j  = evalu( phi, a_j ) #??

    if ( ( phi_j > phi_0 + c1*a_j*dphi_0 ) or (phi_j >= phi_lo) ):
      a_hi = a_j
    else:
      dphi_j = evalu ( dphi, a_j ) #??
      if abs(dphi_j) <= -1*c2*dphi_0 :
        return a_j
      if dphi_j*(a_hi - a_lo) >= 0:
        a_hi = a_lo
      a_lo = a_j
  return

if __name__ == "__main__":
  f,J,H,x1,x2 = rosenbrock()
  X             = np.transpose( np.matrix ([[1.2,1.2]]))
  #X             = np.transpose( np.matrix ([[-1.2,1.]]))
  err           = 0.0001

  x_minimizer, x_points, y_points = grad_desc(f,X,J,H,x1,x2,err)
