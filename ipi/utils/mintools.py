"""Contains simple helper algorithms for minimization.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.

Functions: 
    bracket: Determines the 3 points that bracket the function minimum
    min_brent:  Does one-D minimization (line search) based on bisection method with derivatives. Uses 'bracket' function.
    min_approx: Does approximate n-D minimization (line search) based on sufficient function decrease in the search direction
    BFGS: Constructs an approximate inverse Hessian to determine new search directions. Minimizes the function using 'min_approx' function.
    
"""

__all__ = [ "min_brent" ]

import numpy as np
import math
from ipi.utils.messages import verbosity, warning, info

def bracket(fdf, fdf0=None, x0=0.0, init_step=1.0e-3): 

  """Given an initial point, determines the initial bracket for the minimum
   Arguments:
      fdf = function to minimize [0], derivative of function to minimize [1]
      x0 = initial point
      fdf0 = value of function and its derivative at x0
  """

  # Constants
  gold = 1.618034 # Golden ratio
  glimit = 100.0 # Limit for magnification of parabolic fit step
  tiny = 1.0e-20 # Prevent division by zero

  # x0: initial point of evaluation (e.g. initial atomic position)
  # ax, bx, cx: bracketing points with ax < bx < cx
  # fa, fb, fc: value of function at ax, bx, cx
  if fdf0 is None: fdf0 = fdf(x0)
  ax = x0 
  fa, dfa = fdf0 
  bx = x0 + init_step
  fb, dfb = fdf(bx)
  info(" @BRACKET: Started bracketing", verbosity.debug)
  info(" @BRACKET: Evaluated first step", verbosity.debug)

  # Switch direction to move downhill, if necessary
  if fb > fa:
    tmp = ax
    ax = bx
    bx = tmp
    tmp = fb
    fb = fa
    fa = tmp
    tmp = dfb
    dfb = dfa
    dfa = tmp

  # Initial guess for third bracketing point
  cx = bx + gold * (bx - ax)
  fc, dfc = fdf(cx)
  info(" @BRACKET: Evaluated initial bracket: (%f:%f, %f:%f, %f:%f)" % (ax, fa, bx, fb, cx, fc), verbosity.debug) 

  # Loop until acceptable bracketing condition is achieved
  # u is a point between two of the bracketing points
  # Use parabolic extrapolation to find u. "tiny" prevents possible division by zero
  while fb > fc:
    r = (bx - ax) * (fb - fc)
    q = (bx - cx) * (fb - fa)
    u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * math.copysign(max(abs(q - r), tiny), (q - r))) # Point from parabolic fit
    ulim = bx + glimit * (cx - bx) # Limit for parabolic fit point; *Can test various possibilities*

    # Find minimums between b and c or a and u
    # If parabolic fit unsuccessful, use default step magnification
    # Otherwise:
    # - Parabolic fit between c and its allowed limit
    # - Limit u to maximum allowed value
    # - Use default magnification
    if ((bx - u) * (u - cx)) > 0.0:
      fu, dfu = fdf(u)
      info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
      if fu < fc:
        ax = bx
        bx = u
        fa = fb
        fb = fu
        dfa = dfb
        dfb = dfu
        info(" @BRACKET: Bracketing completed: (%f:%f, %f:%f, %f:%f)" % (ax, fa, bx, fb, cx, fc), verbosity.debug)
        return (ax, bx, cx, fb, dfb)

      elif fu > fb:
        cx = u
        fc = fu
        dfc = dfu
        info(" @BRACKET: Bracketing completed", verbosity.debug)
        return (ax, bx, cx, fb, dfb)

      u = cx + gold * (cx - bx)
      fu, dfu = fdf(u)
      info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
    elif ((cx - u) * (u - ulim)) > 0.0:
      fu, dfu = fdf(u)
      info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
      if fu < fc:
        bx = cx
        cx = u
        u = cx + gold * (cx - bx)
        fb = fc
        fc = fu
        dfb = dfc
        dfc = dfu
        fu, dfu = fdf(u)
        info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
    elif ((u - ulim) * (ulim - cx)) >= 0.0:
      u = ulim
      fu, dfu = fdf(u)
      info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
    else:
      u = cx + gold * (cx - bx)
      fu, dfu = fdf(u)
      info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 

    # Shift points
    ax = bx
    bx = cx
    cx = u
    fa = fb
    fb = fc
    fc = fu
    dfa = dfb
    dfb = dfc
    dfc = dfu

  info(" @BRACKET: Bracketing completed: (%f:%f, %f:%f, %f:%f)" % (ax, fa, bx, fb, cx, fc), verbosity.debug)
  return (ax, bx, cx, fb, dfb)

# One dimensional minimization function using function derivatives
# and Brent's method
def min_brent(fdf, fdf0=None, x0=0.0, tol=1.0e-6, itmax=100, init_step=1.0e-3):

  """Given a maximum number of iterations and a convergence tolerance,
   minimizes the specified function 
   Arguments:
      x0 = initial x-value
      fdf = function to minimize
      fdf0 = initial function value
      tol = convergence tolerance
      itmax = maximum allowed iterations
  """

  # Initializations and constants
  gold = 0.3819660 # Golden ratio
  zeps = 1.0e-10 # Safeguard against trying to find fractional precision for min that is exactly zero
  e = 0.0 # Size of step before last

  # Call initial bracketing routine; takes arguments function and initial x-value
  (ax, bx, cx, fb, dfb) = bracket(fdf, fdf0, x0, init_step)
  
  # Set bracket points
  if ax < cx:
    a = ax
  else:
    a = cx
  if ax > cx:
    b = ax
  else: 
    b = cx

  # Initial points to evaluate
  # g(x) is evaluation of arbitrary function
  # dg(x) is the evaluation of the derivative of g(x)
  x = w = v = bx
  fw = fv = fx = fb # Function
  dfw = dfv = dfx = dfb # Function derivative
  
  # Main loop
  j = 1
  info(" @MINIMIZE: Started 1D minimization", verbosity.debug)
  while j <= itmax:

    # Determine tolerance
    xm = 0.5 * (a + b)
    tol1 = tol * abs(x) + zeps
    tol2 = 2.0 * tol1

    # Test for satisfactory completion
    if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
      info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
      return (x, fx)

    # Initialize d values to outside of bracket
    if abs(e) > tol1:
      d1 = 2.0 * (b - a)
      d2 = d1

      # Secant method with both d points
      if dfw != dfx:
        d1 = (w - x) * dfx / (dfx - dfw)
      if dfv != dfx:
        d2 = (v - x) * dfx / (dfx - dfv)

      # Choose estimate based on derivative at x and move on step
      # before last
      u1 = x + d1
      u2 = x + d2
      ok1 = ((a - u1) * (u1 - b) > 0.0) and (dfx * d1 <= 0.0)
      ok2 = ((a - u2) * (u2 - b) > 0.0) and (dfx * d2 <= 0.0)
      olde = e
      e = d

      # Take an acceptable d; if both are acceptable, choose smallest
      if ok1 or ok2:
        if ok1 and ok2:
          if abs(d1) < abs(d2):
            d = d1
          else:
            d = d2
        elif ok1:
          d = d1
        else:
          d = d2
        if abs (d) <= abs (0.5 * olde):
          u = x + d
          if ((u - a) < tol2) or ((b - u) < tol2):
             d = math.copysign(tol1, (xm - x))
          else:
            if dfx >= 0.0:
              e = a - x
            else:
              e = b - x
            d = 0.5 * e
      else:
        if dfx >= 0.0:
          e = a - x
        else:
          e = b - x
        d = 0.5 * e
    else:
      if dfx >= 0.0:
        e = a - x
      else:
        e = b - x
      d = 0.5 * e
    if abs(d) >= tol1:
      u = x + d
      fu, dfu = fdf(u)
    else:
      u = x + math.copysign(tol1, d)
      fu, dfu = fdf(u)
     
      # If minimum step goes uphill, minimum has been found
      if fu > fx:
        info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
        return (x, fx)
        
    if fu <= fx:
      if u >= x:
        a = x
      else:
        b = x
      v = w
      fv = fw
      dfv = dfw
      w = x
      fw = fx
      dfw = dfx
      x = u
      fx = fu
      dfx = dfu
    else:
      if u < x:
        a = u
      else:
        b = u
      if (fu <= fw) or (w == x):
        v = w
        fv = fw
        dfv = dfw
        w = u
        fw = fu
        dfw = dfu
      elif (fu < fv) or (v == x) or (v == w):
        v = u
        fv = fu
        dfv = dfu
    j += 1
  
  # Exit if maximum number of iterations exceeded
  info(" @MINIMIZE: Error -- maximum iterations for minimization (%d) exceeded, exiting minimization" % itmax, verbosity.low)
  info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
  return (x, fx)

def min_approx(fdf, x0, fdf0=None, d0=None, max_step=100.0, tol=1.0e-6, itmax=100):
  
  """Given an n-dimensional function and its gradient, and an 
  initial point and a direction, finds a new point where the function
  is thought to be 'sufficiently' minimized, i.e. carries out an 
  approximate minimization. 
    Arguments:
      fdf = function and its gradient
      fdf0 = initial function and gradient value
      d0 = n-dimensional initial direction
      p0 = n-dimensional initial point
      max_step = maximum step size
    NOTE: commented loops (for j in range ...) are alternative constructions
    to using NumPy for vector/matrix manipulation
  """
  
  # Initializations and constants
  info(" @MINIMIZE: Started approx. line search for BFGS", verbosity.debug) 
  stepsum = 0.0
  n = len(x0.flatten())
  if fdf0 is None: fdf0 = fdf(x0)
  f0, df0 = fdf0
  if d0 is None: d0 = -df0 / np.sqrt(np.dot(df0.flatten(), df0.flatten()))
  x = np.zeros(n)
  alf = 1.0e-4

  stepsum = np.sqrt(np.dot(d0.flatten(), d0.flatten()))
#  for j in range(0, n):
#    stepsum += d0[j] * d0[j]

  stepsum = np.sqrt(stepsum) 

  # Scale if attempted step is too large
  if stepsum > max_step:
    info(" @MINIMIZE: Scaled step size for line search", verbosity.debug)
    d0 = np.multiply(d0, max_step / stepsum)
#    for j in range(0, n):
#      d0[j] *= max_step / stepsum

#  slope = 0.0

  slope = np.dot(df0.flatten(), d0.flatten())
#  for j in range(0, n):
#    slope += df0[j] * d0[j]

  if slope >= 0.0:
    info(" @MINIMIZE: Warning -- gradient is >= 0 (%f)" % slope, verbosity.low)

  # Compute coefficient for Newton step
#  test = 0.0

  test = np.amax(np.divide(np.absolute(d0), np.maximum(np.absolute(x0), np.ones(n))))
#  for j in range (0, n):
#    tmp = abs(d0[j]) / max(abs(x0[j]), 1.0)
#    if tmp > test:
#      test = tmp

  # Setup to try Newton step first
  alamin = tol / test
  alam = 1.0

  # Minimization Loop
  i = 1
  while i < itmax:
    x = np.add(x0, (alam * d0))
    fx, dfx = fdf(x)
    info(" @MINIMIZE: Calculated energy", verbosity.debug)
#    for j in range(0, n):
#      x[j] = x0[j] + alam * d0[j]
#      fx, dfx = fdf(x)
    
    # Check for convergence on change in x; exit
    if alam < alamin:
      x = x0
      info(" @MINIMIZE: Convergence in position, exited line search", verbosity.debug)
#      for j in range(0, n):
#        x[j] = x0[j]
      return (x, fx)

    # Sufficient function decrease; exit
    elif fx <= (f0 + alf * alam * slope):
      info(" @MINIMIZE: Sufficient function decrease, exited line search", verbosity.debug)
      return (x, fx)

    # No convergence; backtrack
    else:
      info(" @MINIMIZE: No convergence on step; backtrack to find point", verbosity.debug)

      # First backtrack
      if alam == 1.0:
        tmplam = -slope / (2.0 * (fx - f0 - slope))
      
      # Subsequent backtracks
      else:
        rhs1 = fx - f0 - alam * slope
        rhs2 = f2 - f0 - alam2 * slope
        a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2)
        b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2)
        if a == 0.0:
          tmplam = -slope / (2.0 * b)
        
        else:
          disc = b * b - 3.0 * a * slope
          if disc < 0.0:
            tmplam = 0.5 * alam

          elif b <= 0.0:
            tmplam = (-b + np.sqrt(disc)) / (3.0 * a)

          else:
            tmplam = -slope / (b + np.sqrt(disc))

          # Coefficient less than 0.5 * lambda_1
          if tmplam > (0.5 * alam):
            tmplam = 0.5 * alam

    alam2 = alam
    f2 = fx

    # Coefficient greater than 0.1 * lambda_1
    alam = max(tmplam, 0.1 * alam)

    i += 1

  info(" @MINIMIZE: Error - maximum iterations for line search (%d) exceeded, exiting search" % itmax, verbosity.low)
  info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
  return (x, fx)
    
def BFGS(x0, d0, fdf, fdf0=None, invhessian=None, max_step=100, tol=1.0e-6, grad_tol=1.0e-6, itmax=100):
  
  """BFGS minimization. Uses approximate line minimizations.
  Does one step.
    Arguments:
      fdf = function and gradient
      fdf0 = initial function and gradient value
      d0 = initial direction for line minimization
      x0 = initial point
      max_step = limit on step length
      tol = convergence tolerance
      itmax = maximum number of allowed iterations
    NOTE: commented loops (for j in range ...) are alternative constructions
    to using NumPy for vector/matrix manipulation
  """
  
  # Original function value, gradient, other initializations
  zeps = 1.0e-10
  if fdf0 is None: fdf0 = fdf(x0)
  f0, df0 = fdf0
  n = len(x0.flatten())
  if invhessian is None: invhessian = np.eye(n)
  dg = np.zeros(n)
  g = df0.flatten()
  hdg = np.zeros(n)
  x = np.zeros(n)
  linesum = np.dot(x0.flatten(), x0.flatten())
  
  # Initial line direction
  xi = d0

#  for j in range(0, n):
#    for k in range(0, n):
#      invhessian[j][k] = 1.0
#    xi[j] = -g[j]
#    linesum += x0[j] * x0[j]

  # Maximum step size
  max_step = max_step * max(np.sqrt(linesum), n)

 # i = 0 
  #while i <= itmax:

  # Perform approximate line minimization in direction d0
  x, fx = min_approx(fdf, x0, fdf0, xi, max_step, tol, itmax) # must return vectors

  info(" @MINIMIZE: Started BFGS", verbosity.debug)

  # Update line direction (xi) and current point (x0)
  xi = np.subtract(x, x0).flatten()
  x0 = x

    #for j in range(0, n):
     # xi[j] = x[j] - x0[j]
     # x0[j] = x[j]

  # Test for convergence on step size
#  test = 0.0
  test = np.amax(np.divide(np.absolute(xi), np.maximum(np.absolute(x0), np.ones(n))))
#    for j in range(0, n):
#      tmp = abs(xi[j]) / max(abs(x0[j]), 1.0))
#      if tmp > test:
#        test = tmp

  if test < tol:
    info(" @MINIMIZE: Convergence on tolerance, exited BFGS, energy = %f" % fx, verbosity.debug)
    return (x, fx, xi, invhessian)

  # Store old gradient
  dg = g

#    for j in range(0, n):
#      dg[j] = g[j]
      
  # Get new gradient    
  unused, g = fdf(x0)
  info(" @MINIMIZE: Updated gradient", verbosity.debug)
  g = g.flatten()
  test = 0.0
  den = max(fx, 1.0)

  test = np.amax(np.divide(np.multiply(np.absolute(g), np.maximum(np.absolute(x0), np.ones(n))), den))
#    for j in range(0, n):
#      tmp = abs(g[j]) * max(abs(x0[j]), 1.0)) / den
#      if tmp > test:
#        test = tmp

  # Test for convergence on zero gradient
  if test < grad_tol:
    info(" @MINIMIZE: Convergence on zero gradient, exited BFGS, energy = %f" % fx, verbosity.debug)
    return (x, fx, xi, invhessian)

  # Compute difference of gradients
  dg = np.subtract(g, dg)

#    for j in range(0, n):
#      dg[j] = g[j] - dg[j]

  # Difference of gradients times current matrix
#  for j in range(0, n):
#    hdg[j] = 0.0
#    for k in range(0, n):
#      hdg[j] += invhessian[j][k] * dg[k]
  hdg = np.dot(invhessian, dg)
  
#    fac = fae = sumdg = sumxi = 0.0
  fac = np.dot(dg.flatten(), xi.flatten())
  fae = np.dot(dg.flatten(), hdg.flatten())
  sumdg = np.dot(dg.flatten(), dg.flatten())
  sumxi = np.dot(xi.flatten(), xi.flatten())
#    for j in range(0, n):
#      fac += dg[j] * xi[j]
#      fae += dg[j] * hdg[j]
#      sumdg += (dg[j] * dg[j])
#      sumxi += (xi[j] * xi[j])

  # Skip update if not 'fac' sufficiently positive
  if fac > np.sqrt(zeps * sumdg * sumxi):
    info(" @MINIMIZE: Skipped hessian update; direction x gradient insufficient", verbosity.debug)
    fac = 1.0 / fac
    fad = 1.0 / fae

    # Compute BFGS term
    dg = np.subtract(fac * xi, fad * hdg)
#      for j in range(0, n):
#        dg[j] = fac * xi[j] - fad * hdg[j]

#    print "updating hessian"
#    for j in range(0, n):
#      for k in range(j, n):
#        invhessian[j,k] += fac * xi[j] * xi[k] - fad * hdg[j] * hdg[k] + fae * dg[j] * dg[k]
#        invhessian[k,j] = invhessian[j,k]
    invhessian = invhessian + np.outer(xi, xi) * fac - np.outer(hdg, hdg) * fad + np.outer(dg, dg) * fae    
    info(" @MINIMIZE: Updated hessian", verbosity.debug)

  # Update direction
#  for j in range(0, n):
#    xi[j] = 0.0
#    for k in range(0, n):
#      xi[j] -= invhessian[j][k] * g[k]
  xi = np.dot(invhessian, -g)
  info(" @MINIMIZE: Updated search direction", verbosity.debug)
  return (x, fx, xi, invhessian)



