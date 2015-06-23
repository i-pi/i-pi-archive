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
    min_brent:  Applies one-d minimization based on bisection method with derivatives.
    
"""

__all__ = [ "min_brent" ]

import numpy as np
import math
from ipi.utils.messages import verbosity, warning

def bracket(fdf, fdf0=None, x0=0.0, init_step=1e-3): #TODO: ALSO ADD OPTION FOR glimit?

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

  # Loop until acceptable bracketing condition is achieved
  # u is a point between two of the bracketing points
  # Use parabolic extrapolation to find u. "tiny" prevents possible div$
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
      if fu < fc:
        ax = bx
        bx = u
        fa = fb
        fb = fu
        dfa = dfb
        dfb = dfu
        return (ax, bx, cx, fb, dfb)

      elif fu > fb:
        cx = u
        fc = fu
        dfc = dfu
        return (ax, bx, cx, fb, dfb)

      u = cx + gold * (cx - bx)
      fu, dfu = fdf(u)
    elif ((cx - u) * (u - ulim)) > 0.0:
      fu, dfu = fdf(u)
      if fu < fc:
        bx = cx
        cx = u
        u = cx + gold * (cx - bx)
        fb = fc
        fc = fu
        dfb = dfc
        dfc = dfu
        fu, dfu = fdf(u)
    elif ((u - ulim) * (ulim - cx)) >= 0.0:
      u = ulim
      fu, dfu = fdf(u)
    else:
      u = cx + gold * (cx - bx)
      fu, dfu = fdf(u)

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
  return (ax, bx, cx, fb, dfb)

# One dimensional minimization function using function derivatives
# and Brent's method
def min_brent(fdf, fdf0=None, x0=0.0, tol = 1e-6, itmax = 100, init_step = 1e-3):

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
  while j <= itmax:

    # Determine tolerance
    xm = 0.5 * (a + b)
    tol1 = tol * abs(x) + zeps
    tol2 = 2.0 * tol1

    # Test for satisfactory completion
    if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
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
  print "Error: Maximum iterations exceeded:", itmax
  return (x, fx)

