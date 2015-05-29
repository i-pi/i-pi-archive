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


def min_brent(gdg, tol=1e-5, itmax=100, gdg0=None):

  # Initializations and constants
  gold = 0.3819660
  zeps = 1e-10  
  e = 0.0

# MC this is all very confused, you should consider that you are given a starting point
# and find yourself the bracket
  # Requirement that :
  # ax < bx < cx
  # f(bx) < f(ax) and f(bx) < f(cx)
  # near the local minimum
  #!TODO! should find the brackets!
  ax = -1.0
  bx = 0.0
  cx = 1.0

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
  if gdg0 is None: gdg0 = gdg(x)
  fw = fv = fx = gdg0[0]
  dw = dv = dx = gdg0[1]
  
  # Main loop
  iter = 1
  while iter <= itmax:

    # Determine tolerance
    xm = 0.5 * (a + b)
    tol1 = tol * abs(x) + zeps
    tol2 = 2.0 * tol1

    # Test for satisfactory completion
    if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
      xmin = x
      results = [xmin, fx]
      return results

    # Initialize d values to outside of bracket
    if abs(e) > tol1:
      d1 = 2.0 * (b - a)
      d2 = d1

      # Secant method with both d points
      if dw != dx:
        d1 = (w - x) * dx / (dx - dw)
      if dv != dx:
        d2 = (v - x) * dx / (dx - dv)

      # Choose estimate based on derivative at x and move on step
      # before last
      u1 = x + d1
      u2 = x + d2
      ok1 = ((a - u1) * (u1 - b) > 0.0) and (dx * d1 <= 0.0)
      ok2 = ((a - u2) * (u2 - b) > 0.0) and (dx * d2 <= 0.0)
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
             d = abs(tol1) * (xm - x) / abs(xm - x)
          else:
            if dx >= 0.0:
              e = a - x
            else:
              e = b - x
            d = 0.5 * e
      else:
        if dx >= 0.0:
          e = a - x
        else:
          e = b - x
        d = 0.5 * e
    else:
      if dx >= 0.0:
        e = a - x
      else:
        e = b - x
      d = 0.5 * e
    if abs(d) >= tol1:
      u = x + d
      fu, du = gdg(u)
    else:
      u = x + abs(tol1) * d / abs(d)
      fu, du = gdg(u)
     
      # If minimum step goes uphill, minimum has been found
      if fu > fx:
        xmin = x
        return fx
        
    if fu <= fx:
      if u >= x:
        a = x
      else:
        b = x
      v = w
      fv = fw
      dv = dw
      w = x
      fw = fx
      dw = dx
      x = u
      fx = fu
      dx = du
    else:
      if u < x:
        a = u
      else:
        b = u
      if (fu <= fw) or (w == x):
        v = w
        fv = fw
        dv = dw
        w = u
        fw = fu
        dw = du
      elif (fu < fv) or (v == x) or (v == w):
        v = u
        fv = fu
        dv = du
    iter += 1
  
  # Exit if maximum number of iterations exceeded
  print "Error: Maximum iterations exceeded"
  xmin = x
  results = (xmin, fx)
  return results
